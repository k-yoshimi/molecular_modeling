from .build_system_from_toml import build_system_from_toml
from .build_molecule_groups import build_molecule_groups
from .optimize_voronoi_overlap import optimize_voronoi_overlap
from .voronoi_geom import check_system_voronoi_volume_tiling
from .render_voronoi_overlays import render_voronoi_overlays
from .write_structure import write_structure

def main():
    system = build_system_from_toml()
    groups = build_molecule_groups(system)
    system = optimize_voronoi_overlap(system)
    check_system_voronoi_volume_tiling(system, groups, rel_tol=1e-3, abs_tol=1e-6, verbose=True)
    render_voronoi_overlays(system, out_dir="./voronoi_pngs")
    write_structure(system)  
    return

if __name__ == "__main__":
    main()