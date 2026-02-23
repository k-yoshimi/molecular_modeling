from ase.visualize import view
from molecular_modeling.build_system_from_toml import build_system_from_toml
from molecular_modeling.build_molecule_groups import build_molecule_groups
# view(system)
from molecular_modeling_gui import launch_gui
from molecular_modeling.optimize_voronoi_overlap import optimize_voronoi_overlap
from molecular_modeling.voronoi_geom import check_system_voronoi_volume_tiling
from molecular_modeling.render_voronoi_overlays import render_voronoi_overlays
from molecular_modeling.write_structure import write_structure


system = build_system_from_toml(verbose=True)
print("System built")
print("Launching GUI...")
modified_system = launch_gui(system)
print("System modified")
print("Starting optimization...")
groups = build_molecule_groups(system)

rngs = modified_system.info.get("rng", {})
rng_opt = int(rngs.get("optimize", rngs.get("global", 1)))
rng_ren = int(rngs.get("render",   rngs.get("global", 1)))


# Rotate to maximize IoU.
modified_system = optimize_voronoi_overlap(modified_system, groups=groups, rng=rng_opt)

check_system_voronoi_volume_tiling(modified_system, groups, rel_tol=1e-3, abs_tol=1e-6, verbose=True)
render_voronoi_overlays(modified_system, groups=groups, out_dir="./voronoi_pngs", rng=rng_ren)
write_structure(modified_system, outdir="./structure_out", verbose=True)
