import tomllib
import numpy as np
from scipy.spatial import Voronoi
import trimesh
import matplotlib
# Use the 'Agg' backend for headless environments (no display)
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# 1) Load benzene.toml
with open("../../input/benzene.toml", "rb") as f:
    data = tomllib.load(f)

# 2) Build lattice matrix (columns are a, b, c vectors)
lvecs = np.array(data["lattice_vectors"])    # shape=(3,3), rows are vectors
B = lvecs.T                                  # shape=(3,3)
B_inv = np.linalg.inv(B)

# 3) Get center‐of‐mass positions in Cartesian and fractional
coms_cart = np.array(data["center_of_mass"])              # (M,3)
coms_frac = (B_inv @ coms_cart.T).T                        # fractional coords

# 4) Create periodic images: offsets in {-1,0,1}³
offsets = np.array([[i,j,k] for i in (-1,0,1)
                             for j in (-1,0,1)
                             for k in (-1,0,1)])
points_frac = []
center_indices = []
for off_idx, off in enumerate(offsets):
    for pt_idx, frac in enumerate(coms_frac):
        points_frac.append(frac + off)
        # remember index of the central‐cell points (offset == [0,0,0])
        if np.all(off == 0):
            center_indices.append(off_idx * len(coms_frac) + pt_idx)
points_frac = np.vstack(points_frac)                  # (27M,3)
points_cart = (B @ points_frac.T).T                   # (27M,3)

# 5) Compute Voronoi on the expanded point set
vor = Voronoi(points_cart)

# 6) Build unit‐cell box as a convex hull mesh
corners_frac = np.array([[i,j,k] for i in (-1,2)
                                for j in (-1,2)
                                for k in (-1,2)])
corners_cart = (B @ corners_frac.T).T                  # (8,3)
cell_box = trimesh.convex.convex_hull(corners_cart)

# 7) Extract and clip each central Voronoi cell
cells = []
for pt_idx in center_indices:
    region_idx = vor.point_region[pt_idx]
    vert_ids = vor.regions[region_idx]
    if -1 in vert_ids:
        # unbounded region → skip (shouldn't happen with enough padding)
        continue
    verts = vor.vertices[vert_ids]
    # build convex‐hull mesh of the region
    hull = trimesh.convex.convex_hull(verts)
    # intersect with unit cell
    clipped = hull.intersection(cell_box)
    cells.append(clipped)

# 8) Output: volumes of each molecule’s Voronoi cell
for i, cell in enumerate(cells, start=1):
    print(f"Molecule {i}: Voronoi cell volume = {cell.volume:.6f} Å³")

# Optional: visualize the first cell
# cells[0].show()

# Create a new figure with 3D axes
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(projection='3d')

# extract the atom positions list
positions = np.array(data["positions"])

# --- Prepare figure & axes ---
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(projection='3d')

# --- Plot each Voronoi cell ---
for cell in cells:
    verts = cell.vertices                  # (N,3) array
    faces = verts[cell.faces]              # (M,3,3) array
    mesh = Poly3DCollection(
        faces,
        alpha=0.3,                         # transparency
        edgecolor='k',                     # draw edges in black
        linewidth=0.5
    )
    ax.add_collection3d(mesh)

# --- Overlay atomic positions ---
# assume first 6 points are carbons, next 6 are hydrogens
carbons  = positions[:6]
hydrogens = positions[6:]

# plot carbons as larger black spheres
ax.scatter(
    carbons[:,0], carbons[:,1], carbons[:,2],
    s=60,            # marker size
    c='k',           # black
    marker='o',
    label='C'
)
# plot hydrogens as smaller red spheres
ax.scatter(
    hydrogens[:,0], hydrogens[:,1], hydrogens[:,2],
    s=30,            # marker size
    c='r',           # red
    marker='o',
    label='H'
)

# --- Set axis limits to unit cell bounds ---
all_pts = np.vstack([cell.vertices for cell in cells])
xmin, ymin, zmin = all_pts.min(axis=0)
xmax, ymax, zmax = all_pts.max(axis=0)
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_zlim(zmin, zmax)

# --- Labels and legend ---
ax.set_xlabel('X (Å)', fontsize=16, labelpad=10)
ax.set_ylabel('Y (Å)', fontsize=16, labelpad=10)
ax.set_zlabel('Z (Å)', fontsize=16, labelpad=10)

ax.zaxis.set_label_coords(1.05, 0.5)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.legend(prop={'size':18})

plt.subplots_adjust(
    left   = 0.05,   # left margin
    right  = 0.85,   # right margin
    bottom = 0.05,   # bottom margin
    top    = 0.95    # top margin
)

plt.savefig(
    'voronoi_with_atoms.png',
    dpi         = 300,
    bbox_inches = 'tight',
    pad_inches  = 0.6
)

# plt.show()  # if you ever run with a display