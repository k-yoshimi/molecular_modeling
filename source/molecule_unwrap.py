"""Rejoin molecules that are split across periodic boundaries."""

from __future__ import annotations

from collections import deque
from collections.abc import Sequence

import numpy as np
from ase import Atoms
from ase.geometry import find_mic
from ase.neighborlist import natural_cutoffs


def _bond_matrix(
    positions: np.ndarray,
    symbols: Sequence[str],
    cell: np.ndarray | None,
    pbc: bool | Sequence[bool],
) -> tuple[np.ndarray, np.ndarray]:
    """Return (bonded boolean matrix, pair vectors D[i, j] = r_j - r_i with minimum image if periodic)."""
    pos = np.asarray(positions, dtype=float)
    n = len(pos)
    diff = pos[None, :, :] - pos[:, None, :]
    if cell is not None:
        vecs, dist = find_mic(diff.reshape(-1, 3), np.asarray(cell, dtype=float), pbc)
        diff = vecs.reshape(n, n, 3)
        dist = dist.reshape(n, n)
    else:
        dist = np.linalg.norm(diff, axis=2)
    radii = np.array(natural_cutoffs(Atoms(list(symbols)), mult=1.2))
    bonded = dist <= (radii[:, None] + radii[None, :])
    np.fill_diagonal(bonded, False)
    return bonded, diff


def unwrap_molecule(
    positions: np.ndarray,
    symbols: list[str],
    cell: np.ndarray,
    pbc: bool | Sequence[bool] = True,
) -> np.ndarray:
    """Place atoms so that bonded atoms are adjacent, using a breadth-first walk.

    The first atom is kept fixed; every other atom is placed at the minimum-image
    position relative to the already-placed neighbour that reached it. Molecules
    spanning several cell boundaries are handled because displacements accumulate
    along the bond graph.

    Parameters
    ----------
    positions : (N, 3) Cartesian coordinates.
    symbols : chemical symbols.
    cell : (3, 3) cell matrix.
    pbc : periodic flags.

    Returns
    -------
    (N, 3) Cartesian coordinates. Atoms unreachable from atom 0 (disconnected
    fragments) are left at their input positions.
    """
    pos = np.asarray(positions, dtype=float)
    n = len(pos)
    out = pos.copy()
    if n == 0:
        return out
    bonded, diff = _bond_matrix(pos, symbols, np.asarray(cell, dtype=float), pbc)
    seen = np.zeros(n, dtype=bool)
    seen[0] = True
    queue = deque([0])
    while queue:
        i = queue.popleft()
        for j in np.nonzero(bonded[i] & ~seen)[0]:
            out[j] = out[i] + diff[i, j]
            seen[j] = True
            queue.append(j)
    return out


def extract_molecule(atoms: Atoms, indices: Sequence[int]) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Extract a molecule, rejoined across periodic boundaries, centred on its plain-mean centre.

    Returns
    -------
    (symbols, centred positions, centre) where centre is the unweighted mean of the unwrapped positions.
    """
    idx = list(indices)
    symbols = [atoms.get_chemical_symbols()[i] for i in idx]
    pos = atoms.get_positions()[idx]
    unwrapped = unwrap_molecule(pos, symbols, np.array(atoms.get_cell()), atoms.get_pbc())
    com = unwrapped.mean(axis=0)
    return symbols, unwrapped - com, com


def is_connected(positions: np.ndarray, symbols: list[str]) -> bool:
    """Return True if the bond graph (no periodic boundaries) is a single connected component."""
    n = len(positions)
    if n <= 1:
        return True
    bonded, _ = _bond_matrix(np.asarray(positions, dtype=float), symbols, None, False)
    seen = {0}
    stack = [0]
    while stack:
        i = stack.pop()
        for j in np.nonzero(bonded[i])[0]:
            if int(j) not in seen:
                seen.add(int(j))
                stack.append(int(j))
    return len(seen) == n
