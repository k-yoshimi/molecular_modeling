"""Candidate generation from a (substituted) parent crystal."""
from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read
from molecular_modeling.molecule_unwrap import _bond_matrix, is_connected, unwrap_molecule

from .candidate import Candidate
from .config import SeedSpec
from .molecule import match_and_align, molecule_volume


def _check_occupancies(atoms: Atoms, cif: Path) -> None:
    """Raise ValueError if any site of the CIF has an occupancy below 1."""
    occ = atoms.info.get("occupancy")
    if not occ:
        return
    labels = atoms.info.get("_atom_site_label")
    for site, species in occ.items():
        total = sum(species.values())
        if total < 1.0 - 1e-6 or len(species) > 1:
            idx = int(site)
            label = labels[idx] if labels is not None and idx < len(labels) else site
            raise ValueError(
                f"{cif}: site {label} has partial occupancy {dict(species)}; "
                "remove disorder (e.g. in VESTA) and save as P1 before using it as a seed"
            )


def read_seed_molecules(cif: Path) -> tuple[np.ndarray, list[Atoms]]:
    """Read a seed CIF and split it into whole molecules.

    Molecules are the connected components of the bond graph under periodic boundary
    conditions, each rejoined with :func:`unwrap_molecule`.

    Returns
    -------
    (cell, molecules)
        ``cell`` is the (3, 3) cell matrix; each molecule is a non-periodic
        ``ase.Atoms`` in the coordinates of the unwrapped cell.

    Raises
    ------
    ValueError
        If any site has an occupancy below 1 (the message names the site label), or a
        molecule cannot be rejoined into one connected piece.
    """
    cif = Path(cif)
    atoms = read(str(cif), format="cif", store_tags=True)
    _check_occupancies(atoms, cif)  # type: ignore[arg-type]
    cell = np.array(atoms.get_cell())  # type: ignore[union-attr]
    symbols = atoms.get_chemical_symbols()  # type: ignore[union-attr]
    pos = atoms.get_positions()  # type: ignore[union-attr]
    bonded, _ = _bond_matrix(pos, symbols, cell, True)

    seen = np.zeros(len(atoms), dtype=bool)  # type: ignore[arg-type]
    mols: list[Atoms] = []
    for start in range(len(atoms)):  # type: ignore[arg-type]
        if seen[start]:
            continue
        comp = [start]
        seen[start] = True
        stack = [start]
        while stack:
            i = stack.pop()
            for j in np.nonzero(bonded[i] & ~seen)[0]:
                seen[j] = True
                comp.append(int(j))
                stack.append(int(j))
        comp.sort()
        sym = [symbols[i] for i in comp]
        unwrapped = unwrap_molecule(pos[comp], sym, cell, True)
        if not is_connected(unwrapped, sym):
            raise ValueError(f"{cif}: molecule containing atom {comp[0]} cannot be rejoined (bonds wrap around the cell)")
        mols.append(Atoms(sym, positions=unwrapped, pbc=False))
    return cell, mols


def _enumerate_choices(symbols: list[str], substitutions: list[tuple[str, str, int]]) -> list[list[int]]:
    """All selections of atom indices (flattened over rules) in the parent molecule satisfying the rules."""
    per_rule = []
    for old, _new, count in substitutions:
        idx = [i for i, s in enumerate(symbols) if s == old]
        combos = list(itertools.combinations(idx, count))
        if not combos:
            raise ValueError(f"cannot substitute {count} {old} atom(s): the parent molecule has {len(idx)}")
        per_rule.append(combos)
    choices = []
    for pick in itertools.product(*per_rule):
        flat = [i for c in pick for i in c]
        if len(set(flat)) == len(flat):
            choices.append(flat)
    if not choices:
        raise ValueError("no consistent set of substitution choices (rules overlap)")
    return choices


def generate_from_seed(mol: Atoms, seed: SeedSpec) -> list[Candidate]:
    """Overlay ``mol`` rigidly onto every molecule of a substituted parent crystal.

    For each choice of substituted atoms in parent molecule 0, ``mol`` is superposed on the
    substituted parent molecule 0; the rigid transform from parent molecule 0 to each other
    parent molecule is then applied so all molecules are oriented consistently. The cell is
    the parent's, optionally scaled isotropically (``seed.volume_scale``) with molecules
    moved rigidly so their centres keep their fractional coordinates. Atoms are never scaled.

    Returns
    -------
    list of Candidate
        One per substitution choice (one if ``seed.substitutions`` is empty).

    Raises
    ------
    ValueError
        For disordered or fragmented seeds, mismatched parent molecules, or compositions of
        ``mol`` that do not match the substituted parent.
    """
    cell, parents = read_seed_molecules(seed.cif)
    p0 = parents[0]
    p0_sym = p0.get_chemical_symbols()
    for k, p in enumerate(parents[1:], start=1):
        if sorted(p.get_chemical_symbols()) != sorted(p0_sym):
            raise ValueError(f"{seed.cif}: parent molecule {k} has a different composition from molecule 0")

    # parent molecule 0 -> parent molecule k (unsubstituted parents, identical composition)
    k_align = [match_and_align(p0, p, heavy_only=False) for p in parents[1:]]

    if seed.volume_scale == "auto":
        ratio = molecule_volume(mol) / molecule_volume(p0)
    else:
        ratio = float(seed.volume_scale)
    if ratio <= 0:
        raise ValueError(f"volume_scale must be positive, got {seed.volume_scale!r}")
    s = ratio ** (1.0 / 3.0)
    new_cell = cell * s
    inv_cell = np.linalg.inv(cell)

    out: list[Candidate] = []
    for choice in _enumerate_choices(p0_sym, seed.substitutions) if seed.substitutions else [[]]:
        sub_sym = list(p0_sym)
        rule_iter = iter(choice)
        for _old, new, count in seed.substitutions:
            for _ in range(count):
                sub_sym[next(rule_iter)] = new
        p0_sub = Atoms(sub_sym, positions=p0.get_positions(), pbc=False)
        try:
            al = match_and_align(mol, p0_sub, heavy_only=True)
        except ValueError as exc:
            raise ValueError(f"input molecule does not match the substituted parent molecule: {exc}") from exc
        # input molecule in parent-0 atom order, overlaid on parent molecule 0
        base = (mol.get_positions()[al.perm]) @ al.rotation.T + al.translation

        mol_positions = [base]
        for ka in k_align:
            moved = base[ka.perm] @ ka.rotation.T + ka.translation
            mol_positions.append(moved)
        mol_symbols = [sub_sym] + [[sub_sym[j] for j in ka.perm] for ka in k_align]

        all_pos: list[np.ndarray] = []
        all_sym: list[str] = []
        ids: list[int] = []
        for k, (pos, sym) in enumerate(zip(mol_positions, mol_symbols)):
            com = pos.mean(axis=0)
            frac = com @ inv_cell
            pos = pos + (frac @ new_cell - com)
            all_pos.append(pos)
            all_sym.extend(sym)
            ids.extend([k] * len(pos))
        atoms = Atoms(all_sym, positions=np.vstack(all_pos), cell=new_cell, pbc=True)
        atoms.set_array("mol_id", np.array(ids, dtype=int))
        out.append(
            Candidate(
                atoms=atoms,
                origin=f"seed:{Path(seed.cif).stem}",
                space_group=1,
                meta={"substitution_choice": list(choice), "volume_ratio": ratio},
            )
        )
    return out
