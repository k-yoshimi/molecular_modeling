"""Candidate generation from a (substituted) parent crystal."""
from __future__ import annotations

import itertools
import logging
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read
from molecular_modeling.molecule_unwrap import _bond_matrix, is_connected, unwrap_molecule

from .candidate import Candidate
from .config import SeedSpec
from .molecule import AlignResult, match_and_align, molecule_volume

logger = logging.getLogger(__name__)

#: Symmetry detection tolerance (Å) and atom-matching tolerance (Å) for parent-molecule mappings.
SYMPREC = 0.1
MATCH_TOL = 0.3


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
    assert isinstance(atoms, Atoms)
    _check_occupancies(atoms, cif)
    return np.array(atoms.get_cell()), split_molecules(atoms, source=str(cif))


def bond_components(atoms: Atoms) -> list[list[int]]:
    """Atom indices of the connected components of the periodic bond graph of ``atoms``.

    Bonds use the minimum-image convention (same criterion as :func:`unwrap_molecule`). Each
    component is sorted; components are ordered by their first atom.
    """
    bonded, _ = _bond_matrix(atoms.get_positions(), atoms.get_chemical_symbols(), np.array(atoms.get_cell()), True)
    seen = np.zeros(len(atoms), dtype=bool)
    comps: list[list[int]] = []
    for start in range(len(atoms)):
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
        comps.append(sorted(comp))
    return comps


def split_molecules(atoms: Atoms, source: str = "structure") -> list[Atoms]:
    """Split a periodic structure into whole molecules (bond-graph components, unwrapped across PBC).

    Atom order inside each molecule follows the order in ``atoms``; ``mol_id`` is not used.

    Raises
    ------
    ValueError
        If a molecule cannot be rejoined into one connected piece (``source`` names the input).
    """
    cell = np.array(atoms.get_cell())
    symbols = atoms.get_chemical_symbols()
    pos = atoms.get_positions()
    mols: list[Atoms] = []
    for comp in bond_components(atoms):
        sym = [symbols[i] for i in comp]
        unwrapped = unwrap_molecule(pos[comp], sym, cell, True)
        if not is_connected(unwrapped, sym):
            raise ValueError(
                f"{source}: molecule containing atom {comp[0]} cannot be rejoined (bonds wrap around the cell)"
            )
        mols.append(Atoms(sym, positions=unwrapped, pbc=False))
    return mols


def _symmetry_maps(cell: np.ndarray, parents: list[Atoms]) -> list[tuple[np.ndarray, np.ndarray, np.ndarray] | None]:
    """Isometries taking parent molecule 0 onto each parent molecule, from the crystal's own symmetry.

    For each molecule k >= 1 returns ``(M, shift, jmap)`` such that atom ``i`` of molecule 0 at
    Cartesian ``x`` goes to ``x @ M + shift`` and corresponds to atom ``jmap[i]`` of molecule k;
    ``M`` is orthogonal (proper or improper). ``None`` if no space-group operation (plus lattice
    translation) maps molecule 0 onto molecule k.
    """
    import spglib

    syms = [s for p in parents for s in p.get_chemical_symbols()]
    numbers = Atoms(syms).get_atomic_numbers()
    inv = np.linalg.inv(cell)
    frac = [p.get_positions() @ inv for p in parents]
    data = spglib.get_symmetry((cell, np.vstack(frac), numbers),  # type: ignore[arg-type]
         symprec=SYMPREC)
    out: list[tuple[np.ndarray, np.ndarray, np.ndarray] | None] = [None]
    if data is None:
        return out * len(parents)
    rots, trans = data["rotations"], data["translations"]
    z0 = parents[0].get_atomic_numbers()
    for k in range(1, len(parents)):
        zk = parents[k].get_atomic_numbers()
        found = None
        for w_mat, w_vec in zip(rots, trans):
            t = frac[0] @ w_mat.T + w_vec
            n = np.round(frac[k].mean(axis=0) - t.mean(axis=0))
            t = t + n
            cart_t = t @ cell
            cart_k = frac[k] @ cell
            d = np.linalg.norm(cart_t[:, None, :] - cart_k[None, :, :], axis=2)
            d[z0[:, None] != zk[None, :]] = np.inf
            jmap = d.argmin(axis=1)
            if d[np.arange(len(jmap)), jmap].max() < MATCH_TOL and len(set(jmap)) == len(jmap):
                m = inv @ w_mat.T @ cell
                found = (m, (w_vec + n) @ cell, jmap)
                break
        out.append(found)
    return out


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

    # parent molecule 0 -> parent molecule k: crystal symmetry first, rigid fit only for P1 parents
    sym_maps = _symmetry_maps(cell, parents)
    k_align: list[AlignResult | None] = []
    for k in range(1, len(parents)):
        if sym_maps[k] is None:
            logger.warning("%s: no symmetry operation maps molecule 0 onto molecule %d; using a rigid fit", seed.cif, k)
            k_align.append(match_and_align(p0, parents[k], heavy_only=False))
        else:
            k_align.append(None)

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
        mol_symbols = [sub_sym]
        for k, ka in enumerate(k_align, start=1):
            if ka is not None:
                mol_positions.append(base[ka.perm] @ ka.rotation.T + ka.translation)
                mol_symbols.append([sub_sym[j] for j in ka.perm])
            else:
                m, shift, jmap = sym_maps[k]  # type: ignore[misc]
                moved = np.empty_like(base)
                moved[jmap] = base @ m + shift
                syms_k = [""] * len(jmap)
                for i, j in enumerate(jmap):
                    syms_k[j] = sub_sym[i]
                mol_positions.append(moved)
                mol_symbols.append(syms_k)

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
