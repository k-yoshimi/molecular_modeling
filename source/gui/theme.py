# source/gui/theme.py
"""
GUI theme configuration for Molecular Modeling GUI.

Provides centralized color, font, and spacing constants plus a
``configure_theme(root)`` function that applies a ttk theme and
registers custom styles.
"""

import sys
import tkinter as tk
from tkinter import ttk

# ---------------------------------------------------------------------------
# Color palette
# ---------------------------------------------------------------------------

COLORS = {
    # Primary text
    "fg_primary": "#333333",
    "fg_secondary": "#555555",
    "fg_muted": "#888888",
    "fg_hint": "blue",

    # Section headings (collapsible panels)
    "section_fg": "#333333",

    # Active / highlight states
    "bg_active": "#cce5ff",         # drag-mode ON background
    "bg_unfocus": "#ffcccc",        # unfocus button highlight

    # Tooltip
    "tooltip_bg": "#ffffe0",
    "tooltip_fg": "#333333",

    # Trackball canvas
    "tb_bg": "#f0f0f0",
    "tb_border": "#999999",
    "tb_ball_outline": "#888888",
    "tb_ball_fill": "#e8e8e8",
    "tb_inner_fill": "#f4f4f4",
    "tb_crosshair": "#bbbbbb",
    "tb_center_dot": "#999999",
    "tb_active_xy": "#4090e8",
    "tb_active_z": "#e8a040",

    # Background color presets (canvas BG picker)
    "bg_white": "#ffffff",
    "bg_black": "#000000",
    "bg_dark": "#333333",
    "bg_light": "#d0d0d0",
}

# ---------------------------------------------------------------------------
# Font specifications
# ---------------------------------------------------------------------------

FONTS = {
    "section_heading": ("TkDefaultFont", 10, "bold"),
    "body": ("TkDefaultFont", 9),
    "icon_large": ("", 14),
    "icon_medium": ("", 13),
    "icon_small": ("", 8),
    "fixed": ("TkFixedFont", 10),
    "fixed_small": ("Courier", 10),
    "tooltip": ("TkDefaultFont", 9),
}

# ---------------------------------------------------------------------------
# Spacing constants (pixels)
# ---------------------------------------------------------------------------

SPACING = {
    "panel_padx": 8,
    "panel_pady": 4,
    "section_padx": 4,
    "section_pady": 2,
    "button_padx": 4,
    "button_pady": 2,
    "dialog_pad": 10,
    "separator_padx": 4,
    "separator_pady": 2,
    "tooltip_padx": 6,
    "tooltip_pady": 3,
}


def configure_theme(root: tk.Tk) -> None:
    """Apply a ttk theme and register custom styles.

    Parameters
    ----------
    root : tk.Tk
        The root (or toplevel) Tk window.
    """
    style = ttk.Style(root)

    # Pick theme: aqua on macOS (native), clam elsewhere (customisable)
    if sys.platform == "darwin":
        try:
            style.theme_use("aqua")
        except tk.TclError:
            style.theme_use("clam")
    else:
        try:
            style.theme_use("clam")
        except tk.TclError:
            pass  # fallback to default

    # --- Custom label styles ---
    # Do NOT set foreground here — let the system theme handle text color
    # so labels remain readable in both light and dark modes.
    style.configure(
        "SectionHeading.TLabel",
        font=FONTS["section_heading"],
    )
    style.configure(
        "Muted.TLabel",
        font=FONTS["body"],
    )
    style.configure(
        "Hint.TLabel",
        font=FONTS["body"],
    )
