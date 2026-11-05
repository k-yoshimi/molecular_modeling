# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Molecular modeling toolkit built on ASE (Atomic Simulation Environment) with Voronoi cell optimization. Reads TOML-based molecular definitions, builds supercells, optimizes molecular orientations via Voronoi IoU maximization, and exports to QE/VASP/CIF formats. Documentation is primarily in Japanese.

## Commands

All commands run from the repository root with the venv activated:

```bash
# Setup
source .venv/bin/activate
pip install -e ".[dev]"       # install with dev dependencies

# Run the full pipeline
python run.py

# Run tests
pytest                                    # all tests
pytest tests/test_voronoi_geom.py         # single file
pytest tests/test_voronoi_geom.py::TestPolyhedron::test_volume  # single test
pytest -m "not slow"                      # skip slow tests
pytest --cov=molecular_modeling --cov-report=term-missing  # with coverage

# Lint and format
ruff check source/ tests/                # lint
ruff check source/ tests/ --fix          # lint with auto-fix
ruff format source/ tests/               # format
mypy source/ --ignore-missing-imports    # type check (GUI excluded)

# Pre-commit (runs ruff + mypy)
pre-commit run --all-files

# Build docs
SPHINX_LANG=en make -C docs html
SPHINX_LANG=ja make -C docs html

# Launch GUI (ASE-GUI extension with tkinter dialogs)
python molecular_modeling_gui.py

# Crystal structure prediction (needs the csp extra: pip install -e ".[csp]")
molecular-modeling-csp input/csp.toml                     # writes csp_out/
python -m molecular_modeling.csp.validate csp_out exp.cif  # compare with an experimental CIF
```

## Architecture

### Pipeline (`source/main.py`)

The main pipeline runs sequentially:

1. **Build** (`build_system_from_toml.py`) — Parse `input/input.toml` + `data/*.toml`, construct ASE `Atoms` with supercell replication and substitutions. Metadata stored in `system.info` dict.
2. **Group** (`build_molecule_groups.py`) — Index atoms by `mol_id` array into `List[List[int]]` for per-site operations.
3. **Optimize** (`optimize_voronoi_overlap.py`) — Grid search over Euler angles (α,β,γ) to maximize Voronoi cell IoU between system and dictionary cells. Results stored in `system.info` (`optimized_angles_deg`, `optimized_scales`, `optimized_iou`).
4. **Render** (`render_voronoi_overlays.py`) — 3D matplotlib visualization of Voronoi cell overlays → `voronoi_pngs/`.
5. **Export** (`write_structure.py`) — Write `qe.in`, `POSCAR`, `INCAR`, `KPOINTS`, `structure.cif` → `structure_out/`.

### Core Geometry (`voronoi_geom.py`)

The largest module. Contains `Polyhedron` class for convex polyhedra, Voronoi cell computation via scipy halfspace intersection, volume/IoU calculations using trimesh and manifold3d, and periodic boundary wrapping.

### GUI (`source/gui/`)

Tkinter-based dialogs extending ASE-GUI. `controller.py` bridges GUI actions to the pipeline. `molecule_registry.py` manages the molecule database from `data/*.toml`. Excluded from mypy and coverage.

### Crystal Structure Prediction (`source/csp/`)

Separate pipeline that predicts lattice and packing together: generate candidates (`gen_random` PyXtal, `gen_voronoi`, `gen_seed` from a parent crystal) → `screen` → `relax` (GFN2-xTB via tblite, process pool, time budget) → `dedupe` → ranked `candidates.csv` + CIFs. `pipeline.py` orchestrates, `outdir.py` handles staging/lock, `validate.py` compares against experimental CIFs. Tests in `tests/csp/`.

### Input/Output

- **Input config**: `input/input.toml` — cell vectors, `supercell_repeat`, molecule list, `center_of_mass`, substitution rules, RNG seeds
- **Molecule definitions**: `data/<name>.toml` — atomic positions, lattice vectors, center of mass, optional `voronoi_ref_id`
- **Outputs**: `structure_out/` (structure files), `voronoi_pngs/` (visualization)

## Key Conventions

- Python 3.11+ required (configured for 3.13)
- Package mapping: `molecular_modeling` package maps to `source/` directory via `[tool.setuptools] package-dir`
- Import as `from molecular_modeling.<module> import ...` (resolves to `source/<module>.py`)
- pytest configured with `pythonpath = .` in `pytest.ini` so imports resolve from the repository root
- Shared test fixtures in `tests/conftest.py`: `simple_cell`, `simple_atoms`, `simple_polyhedron`
- `system.info` dict on ASE `Atoms` objects carries metadata between pipeline stages
- Ruff rules: E, F, W, I selected; E501 ignored; line length 120
- Coverage excludes `source/gui/` and `source/main.py` (see `.coveragerc`)
