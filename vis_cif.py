#!/usr/bin/env python3
"""
Open a CIF file in the Molecular Modeling GUI.

Usage:
    python vis_cif.py <cif_file>
    python vis_cif.py              # opens a file dialog
"""

import sys

from molecular_modeling_gui import read_cif, launch_gui


def main():
    if len(sys.argv) >= 2:
        atoms = read_cif(sys.argv[1])
    else:
        # No file argument — open a file chooser dialog
        import tkinter as tk
        from tkinter import filedialog

        root = tk.Tk()
        root.withdraw()
        filepath = filedialog.askopenfilename(
            title="Select CIF File",
            filetypes=[("CIF files", "*.cif"), ("All files", "*.*")],
        )
        root.destroy()
        if not filepath:
            print("No file selected.")
            return
        atoms = read_cif(filepath)

    launch_gui(atoms)


if __name__ == "__main__":
    main()
