# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Molecular modeling toolkit built on ASE (Atomic Simulation Environment) with Voronoi cell optimization. Reads TOML-based molecular definitions, builds supercells, optimizes molecular orientations via Voronoi IoU maximization, and exports to QE/VASP/CIF formats. Documentation is primarily in Japanese.

## Commands

All commands run from the repository root with the venv activated:

```bash
# Setup
source .venv/bin/activate

# Run the full pipeline
python run.py

# Run tests
pytest                                    # all tests
pytest tests/test_voronoi_geom.py         # single file
pytest tests/test_voronoi_geom.py::TestPolyhedron::test_volume  # single test
pytest -m "not slow"                      # skip slow tests

# Launch GUI (ASE-GUI extension with tkinter dialogs)
python molecular_modeling_gui.py
```

## Architecture

### Pipeline (`source/main.py`)

The main pipeline runs sequentially:

1. **Build** (`build_system_from_toml.py`) — Parse `input/input.toml` + `data/*.toml`, construct ASE `Atoms` with supercell replication and substitutions. Metadata stored in `system.info` dict.
2. **Group** (`build_molecule_groups.py`) — Index atoms by `mol_id` array into `List[List[int]]` for per-site operations.
3. **Optimize** (`optimize_voronoi_overlap.py`) — Grid search over Euler angles (α,β,γ) to maximize Voronoi cell IoU between system and dictionary cells. Results stored in `system.info` (`optimized_angles_deg`, `optimized_scales`, `optimized_iou`).
4. **Render** (`render_voronoi_overlays.py`) — 3D matplotlib visualization of Voronoi cell overlays → `source/voronoi_pngs/`.
5. **Export** (`write_structure.py`) — Write `qe.in`, `POSCAR`, `INCAR`, `KPOINTS`, `structure.cif` → `source/structure_out/`.

### Core Geometry (`voronoi_geom.py`)

The largest module. Contains `Polyhedron` class for convex polyhedra, Voronoi cell computation via scipy halfspace intersection, volume/IoU calculations using trimesh and manifold3d, and periodic boundary wrapping.

### GUI (`source/gui/`)

Tkinter-based dialogs extending ASE-GUI. `controller.py` bridges GUI actions to the pipeline. `molecule_registry.py` manages the molecule database from `data/*.toml`.

### Input/Output

- **Input config**: `input/input.toml` — cell vectors, `supercell_repeat`, molecule list, `center_of_mass`, substitution rules, RNG seeds
- **Molecule definitions**: `data/<name>.toml` — atomic positions, lattice vectors, center of mass, optional `voronoi_ref_id`
- **Outputs**: `source/structure_out/` (structure files), `source/voronoi_pngs/` (visualization)

## Key Conventions

- Python 3.11+ required (configured for 3.13)
- Pinned dependencies in `pyproject.toml` (setuptools backend)
- Package mapping: `molecular_modeling` package maps to `source/` directory via `[tool.setuptools] package-dir`
- pytest configured with `pythonpath = .` so imports resolve from the repository root
- Shared test fixtures in `tests/conftest.py`: `simple_cell`, `simple_atoms`, `simple_polyhedron`
- `system.info` dict on ASE `Atoms` objects carries metadata between pipeline stages
