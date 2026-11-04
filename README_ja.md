# Molecular Modeling Toolkit

[![Tests](https://github.com/k-yoshimi/molecular_modeling/actions/workflows/test.yml/badge.svg)](https://github.com/k-yoshimi/molecular_modeling/actions/workflows/test.yml)
[![Coverage](https://img.shields.io/endpoint?url=https://gist.githubusercontent.com/k-yoshimi/2905134b2029f53624390ded4fd3ac9c/raw/coverage-badge.json)](https://github.com/k-yoshimi/molecular_modeling)
[![License: LGPL v2.1+](https://img.shields.io/badge/License-LGPL_v2.1+-blue.svg)](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
[![Python >= 3.11](https://img.shields.io/badge/Python-≥3.11-3776AB.svg?logo=python&logoColor=white)](https://www.python.org/)
[![version](https://img.shields.io/badge/version-0.1.0-orange.svg)](https://github.com/k-yoshimi/molecular_modeling)
[![ASE](https://img.shields.io/badge/ASE-3.25-blueviolet.svg)](https://wiki.fysik.dtu.dk/ase/)
[![Doc (en)](https://img.shields.io/badge/docs-English_(PDF)-green)](https://github.com/k-yoshimi/molecular_modeling/blob/main/docs/MolecularModelingToolkit_en.pdf)
[![Doc (ja)](https://img.shields.io/badge/docs-日本語_(PDF)-green)](https://github.com/k-yoshimi/molecular_modeling/blob/main/docs/MolecularModelingToolkit_ja.pdf)

**[English README is here (README.md)](README.md)**

ASE (Atomic Simulation Environment) と Voronoi セル最適化に基づく分子モデリングツールキットです。TOML 形式の分子定義からスーパーセルを構築し、Voronoi IoU 最大化で分子配向を最適化した上で、QE / VASP / CIF 形式の構造ファイルを出力します。

## 主な機能

- **スーパーセル構築**: TOML 設定から ASE `Atoms` オブジェクトを構築、レプリカ・置換に対応
- **Voronoi 最適化**: Euler 角グリッドサーチで Voronoi セルの IoU を最大化し、分子配向を決定
- **構造ファイル出力**: Quantum ESPRESSO (`qe.in`)、VASP (`POSCAR`, `INCAR`, `KPOINTS`)、CIF 形式に対応
- **3D 可視化**: matplotlib による Voronoi セル重ね合わせ画像を生成
- **GUI**: tkinter ベースの対話型分子モデリングインターフェース（ASE-GUI 拡張）
- **結晶構造予測（CSP）**: 孤立分子 1 個から分子結晶の候補を生成（PyXtal によるランダム配置と Voronoi による向き決め、任意で置換元結晶）し、GFN2-xTB で緩和してエネルギー順に並べ、実験 CIF と比較

## セットアップ

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

結晶構造予測を使う場合（PyXtal、tblite、pymatgen が入ります）:

```bash
pip install -e ".[csp]"
```

## 使い方

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

### CSP — 結晶構造予測

```bash
molecular-modeling-csp input/csp.toml                     # 予測（csp_out/ に出力）
python -m molecular_modeling.csp.validate csp_out exp.cif  # 実験 CIF と比較
```

`input/csp.toml` は小さなベンゼンの例です。設定の全項目と出力の読み方は、ドキュメントの CSP ガイドを参照してください。

## ドキュメント

- [日本語ドキュメント (PDF)](docs/MolecularModelingToolkit_ja.pdf)
- [English Documentation (PDF)](docs/MolecularModelingToolkit_en.pdf)

## テスト

```bash
pytest                             # 全テスト
pytest -m "not slow"               # slow マーク以外
pytest tests/test_voronoi_geom.py  # 単一ファイル
```

## リポジトリ構成

```
molecular_modeling/              ← リポジトリルート
├── README.md                    ← 英語版 README
├── README_ja.md                 ← 日本語版 README（このファイル）
├── LICENSE                      ← LGPL-2.1-or-later
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
│   ├── molecule_unwrap.py       ← セル境界で切れた分子をつなぎ直す
│   ├── csp/                     ← 結晶構造予測
│   └── gui/                     ← GUI モジュール
├── input/                       ← 実行時設定 (input.toml, csp.toml)
├── data/                        ← 分子辞書 (*.toml)
├── tests/                       ← テストスイート
├── docs/                        ← ドキュメント (Sphinx / PDF)
├── structure_out/               ← 出力: 構造ファイル
└── voronoi_pngs/                ← 出力: 可視化画像
```

## ライセンス

本プロジェクトは **GNU Lesser General Public License v2.1 以降 (LGPL-2.1-or-later)** の下で公開されています。

詳細は [LICENSE](LICENSE) をご覧ください。
