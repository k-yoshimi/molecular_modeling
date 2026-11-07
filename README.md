# Molecular Modeling Toolkit

[![Tests](https://github.com/k-yoshimi/molecular_modeling/actions/workflows/test.yml/badge.svg)](https://github.com/k-yoshimi/molecular_modeling/actions/workflows/test.yml)
[![Coverage](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/k-yoshimi/2905134b2029f53624390ded4fd3ac9c/raw/coverage-badge.json)](https://github.com/k-yoshimi/molecular_modeling)
[![License: LGPL v2.1+](https://img.shields.io/badge/License-LGPL_v2.1+-blue.svg)](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
[![Python >= 3.11](https://img.shields.io/badge/Python-≥3.11-3776AB.svg?logo=python&logoColor=white)](https://www.python.org/)
[![version](https://img.shields.io/badge/version-0.1.0-orange.svg)](https://github.com/k-yoshimi/molecular_modeling)
[![ASE](https://img.shields.io/badge/ASE-3.25-blueviolet.svg)](https://wiki.fysik.dtu.dk/ase/)
[![Doc (en)](https://img.shields.io/badge/docs-English_(PDF)-green)](https://github.com/k-yoshimi/molecular_modeling/blob/main/docs/MolecularModelingToolkit_en.pdf)
[![Doc (ja)](https://img.shields.io/badge/docs-日本語_(PDF)-green)](https://github.com/k-yoshimi/molecular_modeling/blob/main/docs/MolecularModelingToolkit_ja.pdf)

**[日本語版 README はこちら (README_ja.md)](README_ja.md)**

A molecular modeling toolkit built on ASE (Atomic Simulation Environment) with Voronoi cell optimization. Reads TOML-based molecular definitions, builds supercells, optimizes molecular orientations via Voronoi IoU maximization, and exports structure files in QE / VASP / CIF formats.

## Features

- **Supercell Construction**: Build ASE `Atoms` objects from TOML configuration with replication and substitution support
- **Voronoi Optimization**: Determine molecular orientations by maximizing Voronoi cell IoU via Euler angle grid search
- **Structure File Export**: Quantum ESPRESSO (`qe.in`), VASP (`POSCAR`, `INCAR`, `KPOINTS`), and CIF formats
- **3D Visualization**: Generate Voronoi cell overlay images with matplotlib
- **GUI**: Interactive molecular modeling interface based on tkinter (ASE-GUI extension)
- **Crystal Structure Prediction (CSP)**: Generate candidate molecular crystals from one isolated molecule (PyXtal random and Voronoi-oriented packings, optional seed crystals), relax them with GFN2-xTB (default) or the machine-learned potential OrbMol-v2, rank them by energy and compare them with experimental CIFs

## Setup

Requires Python 3.11 or later.

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -e .
```

For development (testing and documentation):

```bash
pip install -e ".[dev]"
```

For crystal structure prediction (installs PyXtal, tblite and pymatgen):

```bash
pip install -e ".[csp]"
```

## Usage

Run all commands from the repository root.

### CLI — Pipeline Execution

```bash
python run.py
```

Pipeline flow:
1. Read `input/input.toml` + `data/*.toml` and build ASE `Atoms`
2. Optimize molecular orientations by maximizing Voronoi cell IoU
3. Save Voronoi cell overlay images to `voronoi_pngs/`
4. Export structure files to `structure_out/`

### GUI — Interactive Interface

```bash
python molecular_modeling_gui.py
```

### CSP — Crystal Structure Prediction

```bash
molecular-modeling-csp input/csp.toml                     # predict; writes csp_out/
python -m molecular_modeling.csp.validate csp_out exp.cif  # compare with an experimental CIF
```

`input/csp.toml` is a small benzene example. See the CSP guide in the documentation for every setting and how to read the output.

## Documentation

- [English Documentation (PDF)](docs/MolecularModelingToolkit_en.pdf)
- [日本語ドキュメント (PDF)](docs/MolecularModelingToolkit_ja.pdf)

## Tests

```bash
pytest                             # all tests
pytest -m "not slow"               # skip slow tests
pytest tests/test_voronoi_geom.py  # single file
```

## Repository Structure

```
molecular_modeling/              <- repository root
├── README.md                    <- this file (English)
├── README_ja.md                 <- Japanese README
├── LICENSE                      <- LGPL-2.1-or-later
├── .github/workflows/           <- CI (GitHub Actions)
├── pyproject.toml               <- dependencies and build config
├── run.py                       <- CLI entry point
├── molecular_modeling_gui.py    <- GUI entry point
├── source/                      <- main package (= molecular_modeling)
│   ├── main.py                  <- main pipeline
│   ├── voronoi_geom.py          <- Voronoi geometry
│   ├── build_system_from_toml.py
│   ├── optimize_voronoi_overlap.py
│   ├── render_voronoi_overlays.py
│   ├── write_structure.py
│   ├── molecule_unwrap.py       <- rejoin molecules split across the cell boundary
│   ├── csp/                     <- crystal structure prediction
│   └── gui/                     <- GUI modules
├── input/                       <- runtime config (input.toml, csp.toml)
├── data/                        <- molecule definitions (*.toml)
├── tests/                       <- test suite
├── docs/                        <- documentation (Sphinx / PDF)
├── structure_out/               <- output: structure files
└── voronoi_pngs/                <- output: visualization images
```

## License

This project is licensed under the **GNU Lesser General Public License v2.1 or later (LGPL-2.1-or-later)**.

See [LICENSE](LICENSE) for details.
