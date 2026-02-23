# celltools_gui_326.py
# ASE 3.26.0 / ASE-GUI 拡張：CellTools メニューに triclinic 対応を追加
# - Auto from extents (orthorhombic + padding)
# - Pad along lattice (keep angles)           ← 角度を保持して a/b/c 方向に余白追加
# - Set triclinic from (a,b,c,α,β,γ)          ← 任意の triclinic セルを指定して反映

import numpy as np
from ase.io import read
from ase.gui.gui import GUI
import ase.gui.ui as ui

# 国際化文字列 _ を ASE から import（フォールバックも用意）
try:
    from ase.gui.i18n import _
except Exception:
    def _(s): return s

# 追加: 入力ダイアログ用（軽量に tkinter の simpledialog を利用）
try:
    import tkinter as tk
    from tkinter import simpledialog, messagebox
except Exception:
    tk = None
    simpledialog = None
    messagebox = None


class CellToolsGUI(GUI):
    # ===== 既存：直交セル + 余白 =====
    def auto_cell_from_extents(self, key=None, pad=(2.0, 2.0, 15.0),
                               pbc=(True, True, False), scale=False):
        """表示中の Atoms から外接直方体＋余白で新しい直交セルを作る"""
        atoms = self.atoms  # 現在のフレーム
        pos = atoms.get_positions()
        mins = pos.min(axis=0)
        maxs = pos.max(axis=0)
        lengths = (maxs - mins) + np.array(pad, dtype=float)
        newcell = np.diag(lengths)  # 直交セル

        atoms.set_pbc(pbc)
        atoms.set_cell(newcell, scale_atoms=scale)

        self.set_frame()
        self.draw()

    # ===== 追加１：角度を保持したまま a/b/c 方向に余白 =====
    def pad_along_lattice_keep_angles(self, key=None,
                                      pad=(0.0, 0.0, 15.0),
                                      pbc=(True, True, False),
                                      scale=False,
                                      ask=True):
        """
        いまのセル角度(α,β,γ)を保持しつつ、各格子方向にパディング[Å]を追加。
        pad は (pad_a, pad_b, pad_c) を「両側」ぶんとして扱うので、長さは +2*pad_*。
        ask=True のときはダイアログで pad を指定可能。
        """
        atoms = self.atoms
        cell = atoms.get_cell()            # 3x3 (a,b,c ベクトル)
        a, b, c = cell[0], cell[1], cell[2]

        if ask and simpledialog is not None:
            root = self.window.win if hasattr(self, "window") else None
            ans = simpledialog.askstring(
                parent=root,
                title="Pad along lattice (Å)",
                prompt="Enter pad_a, pad_b, pad_c in Å (comma-separated):\n例: 0, 0, 15"
            )
            try:
                if ans:
                    pad = tuple(float(x.strip()) for x in ans.split(","))
                    if len(pad) != 3:
                        raise ValueError
            except Exception:
                if messagebox:
                    messagebox.showerror("Error", "Invalid input. Using defaults.")

        # 単位ベクトル
        def _unit(v):
            n = np.linalg.norm(v)
            return v / n if n > 1e-12 else v

        ahat = _unit(a)
        bhat = _unit(b)
        chat = _unit(c)

        new_a = a + 2.0 * pad[0] * ahat
        new_b = b + 2.0 * pad[1] * bhat
        new_c = c + 2.0 * pad[2] * chat

        atoms.set_pbc(pbc)
        atoms.set_cell([new_a, new_b, new_c], scale_atoms=scale)
        self.draw()

    # ===== 追加２：任意の triclinic パラメータを指定 =====
    def set_triclinic_from_cellpar(self, key=None,
                                   a=10.0, b=12.0, c=14.0,
                                   alpha=92.0, beta=104.0, gamma=88.0,
                                   pbc=(True, True, True),
                                   scale=False,
                                   ask=True):
        """
        指定した (a,b,c,α,β,γ) から triclinic セルを生成して反映。
        角度は度（degrees）。ask=True のときはダイアログで数値指定。
        """
        from ase.geometry import cellpar_to_cell, cell_to_cellpar

        atoms = self.atoms

        if ask and simpledialog is not None:
            # 現在セルのパラメータを初期値にして入力を促す
            cur = cell_to_cellpar(atoms.get_cell())
            hint = f"{cur[0]:.3f}, {cur[1]:.3f}, {cur[2]:.3f}, {cur[3]:.3f}, {cur[4]:.3f}, {cur[5]:.3f}"
            root = self.window.win if hasattr(self, "window") else None
            ans = simpledialog.askstring(
                parent=root,
                title="Set triclinic from (a,b,c,α,β,γ)",
                prompt=("Enter a,b,c,alpha,beta,gamma (degrees), comma-separated:\n"
                        f"例: {hint}")
            )
            try:
                if ans:
                    vals = tuple(float(x.strip()) for x in ans.split(","))
                    if len(vals) != 6:
                        raise ValueError
                    a, b, c, alpha, beta, gamma = vals
            except Exception:
                if messagebox:
                    messagebox.showerror("Error", "Invalid input. Using defaults.")

        newcell = cellpar_to_cell([a, b, c, alpha, beta, gamma])  # 3x3
        atoms.set_pbc(pbc)
        atoms.set_cell(newcell, scale_atoms=scale)
        self.draw()

    # ===== メニューを拡張 =====
    def get_menu_data(self):
        """既存メニューに 'CellTools' を追加し、項目を並べる"""
        menu = super().get_menu_data()
        M = ui.MenuItem
        menu.append((
            _('CellTools'),
            [
                # 既存の orthorhombic + padding
                M(_('Auto from extents (orthorhombic + padding)'),
                  self.auto_cell_from_extents),

                # ここから下に追記：triclinic 系の 2 機能
                M(_('Pad along lattice (keep angles)'),
                  self.pad_along_lattice_keep_angles),
                M(_('Set triclinic from (a,b,c,α,β,γ)'),
                  self.set_triclinic_from_cellpar),
            ]
        ))
        return menu


if __name__ == "__main__":
    # 使い方:
    # python celltools_gui_326.py your_structure.cif
    import sys
    if len(sys.argv) < 2:
        print("Usage: python celltools_gui_326.py <structure file (cif/xyz/vasp...)>")
        sys.exit(0)
    atoms = read(sys.argv[1])  # 任意の入力ファイル
    gui = CellToolsGUI(images=[atoms])
    gui.run()
