# Molecular Modeling Toolkit

[![Tests](https://github.com/k-yoshimi/molecular_modeling/actions/workflows/test.yml/badge.svg)](https://github.com/k-yoshimi/molecular_modeling/actions/workflows/test.yml)

ASE (Atomic Simulation Environment) と Voronoi セル最適化に基づく分子モデリングツールキットです。TOML 形式の分子定義からスーパーセルを構築し、Voronoi IoU 最大化で分子配向を最適化した上で、QE / VASP / CIF 形式の構造ファイルを出力します。

A molecular modeling toolkit built on ASE with Voronoi cell optimization. Reads TOML-based molecular definitions, builds supercells, optimizes molecular orientations, and exports to QE/VASP/CIF formats.

## 主な機能 / Features

- **スーパーセル構築**: TOML 設定から ASE `Atoms` オブジェクトを構築、レプリカ・置換に対応
- **Voronoi 最適化**: Euler 角グリッドサーチで Voronoi セルの IoU を最大化し、分子配向を決定
- **構造ファイル出力**: Quantum ESPRESSO (`qe.in`)、VASP (`POSCAR`, `INCAR`, `KPOINTS`)、CIF 形式に対応
- **3D 可視化**: matplotlib による Voronoi セル重ね合わせ画像を生成
- **GUI**: tkinter ベースの対話型分子モデリングインターフェース（ASE-GUI 拡張）

## セットアップ / Setup

Python 3.11 以降が必要です。

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -e .
```

開発用（テスト・ドキュメント生成）:

```bash
pip install -e ".[dev]"
```

## 使い方 / Usage

リポジトリルートで実行してください。

### CLI — パイプライン実行

```bash
python run.py
```

処理フロー:
1. `input/input.toml` + `data/*.toml` を読み込み → ASE `Atoms` 構築
2. Voronoi セル IoU 最大化で分子配向を最適化
3. Voronoi セルの重ね合わせ画像を `voronoi_pngs/` に保存
4. 構造ファイルを `structure_out/` に出力

### GUI — 対話型インターフェース

```bash
python molecular_modeling_gui.py
```

## ドキュメント / Documentation

- [日本語ドキュメント (PDF)](docs/MolecularModelingToolkit_ja.pdf)
- [English Documentation (PDF)](docs/MolecularModelingToolkit_en.pdf)

## テスト / Tests

```bash
pytest                         # 全テスト
pytest -m "not slow"           # slow マーク以外
pytest tests/test_voronoi_geom.py  # 単一ファイル
```

## リポジトリ構成 / Repository Structure

```
molecular_modeling/              ← リポジトリルート
├── README.md
├── CLAUDE.md
├── .github/workflows/           ← CI (GitHub Actions)
├── pyproject.toml               ← 依存関係・ビルド設定
├── run.py                       ← CLI エントリーポイント
├── molecular_modeling_gui.py    ← GUI エントリーポイント
├── source/                      ← パッケージ本体 (= molecular_modeling)
│   ├── main.py                  ← メインパイプライン
│   ├── voronoi_geom.py          ← Voronoi 幾何計算
│   ├── build_system_from_toml.py
│   ├── optimize_voronoi_overlap.py
│   ├── render_voronoi_overlays.py
│   ├── write_structure.py
│   └── gui/                     ← GUI モジュール
├── input/                       ← 実行時設定 (input.toml)
├── data/                        ← 分子辞書 (*.toml)
├── tests/                       ← テストスイート
├── docs/                        ← ドキュメント (Sphinx / PDF)
├── structure_out/               ← 出力: 構造ファイル
└── voronoi_pngs/                ← 出力: 可視化画像
```

## ライセンス / License

This project is licensed under the **GNU Lesser General Public License v2.1 or later (LGPL-2.1-or-later)**.

本プロジェクトは **GNU Lesser General Public License v2.1 以降 (LGPL-2.1-or-later)** の下で公開されています。

詳細は [LICENSE](LICENSE) をご覧ください。See [LICENSE](LICENSE) for details.
