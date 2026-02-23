from .build_molecule_groups import build_molecule_groups
from .build_system_from_toml import build_system_from_toml
from .optimize_voronoi_overlap import optimize_voronoi_overlap
from .render_voronoi_overlays import render_voronoi_overlays
from .voronoi_geom import check_system_voronoi_volume_tiling
from .write_structure import write_structure


def main():
    """
    Entry point: build system, optimize orientations, render, and export structures.

    Pipeline
    --------
    1. Build `Atoms` from TOML inputs.
    2. Group atoms by molecule/site.
    3. Optimize Voronoi IoU via rotations (per site).
    4. Check Voronoi volume tiling consistency.
    5. Render overlays to PNGs.
    6. Write QE/VASP/CIF outputs.
    """
    system = build_system_from_toml(verbose=True)
    groups = build_molecule_groups(system)

    rngs = system.info.get("rng", {})
    rng_opt = int(rngs.get("optimize", rngs.get("global", 1)))
    rng_ren = int(rngs.get("render",   rngs.get("global", 1)))

    # Rotate to maximize IoU.
    system = optimize_voronoi_overlap(system, groups=groups, rng=rng_opt)

    check_system_voronoi_volume_tiling(system, groups, rel_tol=1e-3, abs_tol=1e-6, verbose=True)
    render_voronoi_overlays(system, groups=groups, out_dir="./voronoi_pngs", rng=rng_ren)
    write_structure(system, outdir="./structure_out", verbose=True)
    return

if __name__ == "__main__":
    main()
