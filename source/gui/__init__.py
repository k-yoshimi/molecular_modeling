# source/gui/__init__.py
"""
GUI module for molecular_modeling.

Provides tkinter-based dialogs for ASE-GUI integration.
"""

from .add_molecule_dialog import AddMoleculeDialog, show_add_molecule_dialog
from .base_dialog import BaseDialog, CoordinateEntry, LabeledCombobox, LabeledEntry
from .controller import MolecularController, get_controller, reset_controller
from .export_dialog import ExportDialog, show_export_dialog
from .molecule_list_dialog import MoleculeListDialog, show_molecule_list_dialog
from .molecule_registry import MoleculeRegistry, get_registry
from .substitution_dialogs import (
    AtomSubstitutionDialog,
    GlobalAtomSubstitutionDialog,
    MoleculeAtomSubstitutionDialog,
    MoleculeSubstitutionDialog,
    show_atom_substitution_dialog,
    show_global_substitution_dialog,
    show_molecule_atom_substitution_dialog,
    show_molecule_substitution_dialog,
)
from .supercell_dialog import SupercellDialog, show_supercell_dialog

__all__ = [
    # Base
    "BaseDialog",
    "LabeledEntry",
    "LabeledCombobox",
    "CoordinateEntry",
    # Registry
    "MoleculeRegistry",
    "get_registry",
    # Dialogs
    "AddMoleculeDialog",
    "show_add_molecule_dialog",
    "MoleculeListDialog",
    "show_molecule_list_dialog",
    "GlobalAtomSubstitutionDialog",
    "MoleculeSubstitutionDialog",
    "MoleculeAtomSubstitutionDialog",
    "AtomSubstitutionDialog",
    "show_global_substitution_dialog",
    "show_molecule_substitution_dialog",
    "show_molecule_atom_substitution_dialog",
    "show_atom_substitution_dialog",
    "SupercellDialog",
    "show_supercell_dialog",
    "ExportDialog",
    "show_export_dialog",
    # Controller
    "MolecularController",
    "get_controller",
    "reset_controller",
]
