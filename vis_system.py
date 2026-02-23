from ase.visualize import view
from source.build_system_from_toml import build_system_from_toml
system = build_system_from_toml(verbose=True)
# view(system)
from molecular_modeling_gui import launch_gui
launch_gui(system)
