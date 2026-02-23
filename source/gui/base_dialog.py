# source/gui/base_dialog.py
"""
Base dialog class for molecular modeling GUI dialogs.

Provides common functionality for tkinter dialogs used in ASE-GUI extensions.
"""

import tkinter as tk
from tkinter import messagebox, ttk
from typing import Any, Callable, Optional


class BaseDialog(tk.Toplevel):
    """
    Base class for modal dialogs in molecular modeling GUI.

    Provides common patterns:
    - Modal behavior with grab_set
    - Standard OK/Cancel buttons
    - Result handling
    - Error display via messagebox

    Parameters
    ----------
    parent : tk.Widget
        Parent widget (usually the main GUI window).
    title : str
        Dialog window title.
    """

    def __init__(self, parent: Optional[tk.Widget], title: str = "Dialog"):
        super().__init__(parent)
        self.title(title)
        self.parent = parent
        self.result: Any = None
        self._cancelled = True

        # Keep dialog on top of parent but allow interaction with parent window
        self.transient(parent)

        # Create main frame
        self.main_frame = ttk.Frame(self, padding="10")
        self.main_frame.pack(fill=tk.BOTH, expand=True)

        # Build dialog body (override in subclasses)
        self.body_frame = ttk.Frame(self.main_frame)
        self.body_frame.pack(fill=tk.BOTH, expand=True, pady=(0, 10))
        self.build_body(self.body_frame)

        # Button frame
        self.button_frame = ttk.Frame(self.main_frame)
        self.button_frame.pack(fill=tk.X)
        self.build_buttons(self.button_frame)

        # Bind events
        self.protocol("WM_DELETE_WINDOW", self.on_cancel)
        self.bind("<Escape>", lambda e: self.on_cancel())
        self.bind("<Return>", lambda e: self.on_ok())

        # Center on parent
        self.center_on_parent()

        # Focus
        self.focus_set()

    def build_body(self, parent: ttk.Frame) -> None:
        """
        Build the dialog body. Override in subclasses.

        Parameters
        ----------
        parent : ttk.Frame
            Frame to build body widgets in.
        """
        pass

    def build_buttons(self, parent: ttk.Frame) -> None:
        """
        Build standard OK/Cancel buttons.

        Parameters
        ----------
        parent : ttk.Frame
            Frame to build buttons in.
        """
        cancel_btn = ttk.Button(parent, text="Cancel", command=self.on_cancel)
        cancel_btn.pack(side=tk.RIGHT, padx=(5, 0))

        ok_btn = ttk.Button(parent, text="OK", command=self.on_ok)
        ok_btn.pack(side=tk.RIGHT)

    def validate(self) -> bool:
        """
        Validate dialog input. Override in subclasses.

        Returns
        -------
        bool
            True if input is valid, False otherwise.
        """
        return True

    def apply(self) -> None:
        """
        Apply dialog result. Override in subclasses.

        Sets self.result to the dialog output.
        """
        pass

    def on_ok(self) -> None:
        """Handle OK button click."""
        if not self.validate():
            return
        self._cancelled = False
        self.apply()
        self.destroy()

    def on_cancel(self) -> None:
        """Handle Cancel button click."""
        self._cancelled = True
        self.result = None
        self.destroy()

    def center_on_parent(self) -> None:
        """Center the dialog on its parent window."""
        self.update_idletasks()
        if self.parent:
            parent_x = self.parent.winfo_rootx()
            parent_y = self.parent.winfo_rooty()
            parent_w = self.parent.winfo_width()
            parent_h = self.parent.winfo_height()
            w = self.winfo_width()
            h = self.winfo_height()
            x = parent_x + (parent_w - w) // 2
            y = parent_y + (parent_h - h) // 2
            self.geometry(f"+{x}+{y}")

    @property
    def cancelled(self) -> bool:
        """Return True if dialog was cancelled."""
        return self._cancelled

    @staticmethod
    def show_error(title: str, message: str) -> None:
        """
        Show an error message box.

        Parameters
        ----------
        title : str
            Message box title.
        message : str
            Error message to display.
        """
        messagebox.showerror(title, message)

    @staticmethod
    def show_info(title: str, message: str) -> None:
        """
        Show an info message box.

        Parameters
        ----------
        title : str
            Message box title.
        message : str
            Info message to display.
        """
        messagebox.showinfo(title, message)

    @staticmethod
    def ask_yes_no(title: str, message: str) -> bool:
        """
        Show a yes/no confirmation dialog.

        Parameters
        ----------
        title : str
            Message box title.
        message : str
            Question to display.

        Returns
        -------
        bool
            True if user clicked Yes, False otherwise.
        """
        return messagebox.askyesno(title, message)


class LabeledEntry(ttk.Frame):
    """
    A labeled entry widget for forms.

    Parameters
    ----------
    parent : tk.Widget
        Parent widget.
    label : str
        Label text.
    default : str, optional
        Default entry value.
    width : int, optional
        Entry width in characters.
    """

    def __init__(
        self,
        parent: tk.Widget,
        label: str,
        default: str = "",
        width: int = 20,
    ):
        super().__init__(parent)
        self.label = ttk.Label(self, text=label, width=15, anchor=tk.W)
        self.label.pack(side=tk.LEFT)

        self.var = tk.StringVar(value=default)
        self.entry = ttk.Entry(self, textvariable=self.var, width=width)
        self.entry.pack(side=tk.LEFT, fill=tk.X, expand=True)

    def get(self) -> str:
        """Get the entry value."""
        return self.var.get()

    def set(self, value: str) -> None:
        """Set the entry value."""
        self.var.set(value)


class LabeledCombobox(ttk.Frame):
    """
    A labeled combobox widget for forms.

    Parameters
    ----------
    parent : tk.Widget
        Parent widget.
    label : str
        Label text.
    values : list[str]
        Combobox options.
    default : str, optional
        Default selected value.
    width : int, optional
        Combobox width in characters.
    """

    def __init__(
        self,
        parent: tk.Widget,
        label: str,
        values: list,
        default: str = "",
        width: int = 18,
    ):
        super().__init__(parent)
        self.label = ttk.Label(self, text=label, width=15, anchor=tk.W)
        self.label.pack(side=tk.LEFT)

        self.var = tk.StringVar(value=default if default else (values[0] if values else ""))
        self.combobox = ttk.Combobox(
            self, textvariable=self.var, values=values, width=width, state="readonly"
        )
        self.combobox.pack(side=tk.LEFT, fill=tk.X, expand=True)

    def get(self) -> str:
        """Get the selected value."""
        return self.var.get()

    def set(self, value: str) -> None:
        """Set the selected value."""
        self.var.set(value)

    def set_values(self, values: list) -> None:
        """Update the combobox options."""
        self.combobox["values"] = values
        if values and not self.var.get():
            self.var.set(values[0])


class CoordinateEntry(ttk.Frame):
    """
    A 3D coordinate entry widget (X, Y, Z).

    Parameters
    ----------
    parent : tk.Widget
        Parent widget.
    label : str
        Label text.
    defaults : tuple[float, float, float], optional
        Default X, Y, Z values.
    """

    def __init__(
        self,
        parent: tk.Widget,
        label: str,
        defaults: tuple = (0.0, 0.0, 0.0),
    ):
        super().__init__(parent)
        self.label = ttk.Label(self, text=label, width=15, anchor=tk.W)
        self.label.pack(side=tk.LEFT)

        self.x_var = tk.StringVar(value=str(defaults[0]))
        self.y_var = tk.StringVar(value=str(defaults[1]))
        self.z_var = tk.StringVar(value=str(defaults[2]))

        coord_frame = ttk.Frame(self)
        coord_frame.pack(side=tk.LEFT, fill=tk.X, expand=True)

        ttk.Label(coord_frame, text="X:").pack(side=tk.LEFT)
        ttk.Entry(coord_frame, textvariable=self.x_var, width=8).pack(side=tk.LEFT, padx=(0, 5))

        ttk.Label(coord_frame, text="Y:").pack(side=tk.LEFT)
        ttk.Entry(coord_frame, textvariable=self.y_var, width=8).pack(side=tk.LEFT, padx=(0, 5))

        ttk.Label(coord_frame, text="Z:").pack(side=tk.LEFT)
        ttk.Entry(coord_frame, textvariable=self.z_var, width=8).pack(side=tk.LEFT)

    def get(self) -> tuple:
        """
        Get the coordinate values.

        Returns
        -------
        tuple[float, float, float]
            X, Y, Z values.

        Raises
        ------
        ValueError
            If any value is not a valid float.
        """
        return (
            float(self.x_var.get()),
            float(self.y_var.get()),
            float(self.z_var.get()),
        )

    def set(self, x: float, y: float, z: float) -> None:
        """Set the coordinate values."""
        self.x_var.set(str(x))
        self.y_var.set(str(y))
        self.z_var.set(str(z))
