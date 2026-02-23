# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [0.1.0] - 2025-07-28

### Added
- Core pipeline: build system from TOML, optimize Voronoi overlap, render overlays, export structures
- `build_system_from_toml`: parse `input/input.toml` + `data/*.toml`, construct ASE `Atoms` with supercell replication and substitutions
- `optimize_voronoi_overlap`: grid search + Nelder-Mead over Euler angles to maximize Voronoi IoU
- `render_voronoi_overlays`: 3D matplotlib visualization of Voronoi cell overlays
- `write_structure`: export to QE (`qe.in`), VASP (`POSCAR`, `INCAR`, `KPOINTS`), and CIF formats
- `voronoi_geom`: Voronoi cell computation via scipy halfspace intersection, volume/IoU calculations
- `check_system_voronoi_volume_tiling`: sanity checker for Voronoi cell volume tiling
- Configurable `data_dir` key in `input.toml` for molecule TOML file location
- GUI: tkinter-based molecular modeling interface extending ASE-GUI
- CLI entry point: `molecular-modeling` command
- Comprehensive test suite (261 tests, 98% core logic coverage)
- CI via GitHub Actions (Python 3.11 and 3.13) with coverage badge
- Documentation in English and Japanese (Sphinx PDF)
- LGPL-2.1-or-later license

### Performance
- Cached system seeds/volume in IoU computation (avoids redundant Voronoi calculations during optimization)
- Vectorized halfspace construction, lattice translations, and unit cell box mesh generation
