# molecular_modeling

分子モデリング（ASE + Voronoi 最適化）を Python で実行するためのツール一式です。`molecular_modeling` ディレクトリ内で完結しており、TOML で定義した分子配置を最適化したうえで構造ファイル・可視化画像を生成します。

## セットアップ
- 前提: Python 3.11 以降（`python -V` で確認）。
- 依存関係のインストール:
  ```bash
  cd molecular_modeling
  python -m venv .venv
  source .venv/bin/activate
  pip install --upgrade pip
  pip install -r requirements.txt
  ```
- 主な依存モジュール（`requirements.txt` に全記載）:
  - 計算: `ase`, `numpy`, `scipy`
  - 描画/可視化: `matplotlib`, `contourpy`, `pillow`
  - 形状処理: `trimesh`, `manifold3d`
  - 補助: `packaging`, `python-dateutil`, `cycler`, `kiwisolver`, `fonttools`

## 実行手順
`molecular_modeling` ディレクトリ直下で以下を実行します。
```bash
python run.py
```
主なフローは次の通りです。
1. `input/input.toml` と `data/*.toml` を読み込み、ASE の `Atoms` へ組み立て。
2. Voronoi 体積の重なりを最小化するように分子配向を最適化。
3. Voronoi セルの重ね合わせ画像を保存。
4. 各種計算コード向け入力（QE, VASP, CIF など）を出力。

## 入出力のポイント
- 入力設定: `input/input.toml`
  - `cell` と `supercell_repeat` でスーパーセルを定義。
  - `molecules` に使用する分子辞書名（`data/<name>.toml`）を列挙。
  - `center_of_mass` で各分子のサイト位置を指定（`use_fractional=true` で分率座標解釈）。
  - 置換指定（`molecule_substitutions` / `atom_substitutions` / `global_atom_substitutions` など）で分子・原子種を差し替え可能。
  - 乱数設定: `rng`（または `optimize_rng` / `render_rng`）で最適化・描画の乱数シードを制御。
- 分子辞書: `data/*.toml`
  - `positions` と `center_of_mass`（必要なら `use_fractional=true`）で分子の基準配置を定義。
  - `voronoi_ref_id` を指定すると、その分子を Voronoi 参照として利用します。
- 出力:
  - 構造ファイル: `structure_out/` に `qe.in`, `POSCAR`, `KPOINTS`, `INCAR`, `structure.cif` などを生成。
  - 可視化: `voronoi_pngs/` に Voronoi セルの重ね合わせ画像を保存。

## リポジトリ構成 (抜粋)
- `molecular_modeling/run.py`: 実行エントリーポイント。
- `molecular_modeling/source/`: 実装（セル構築・最適化・描画・出力）。
- `molecular_modeling/input/input.toml`: 実行時設定。
- `molecular_modeling/data/`: 分子ごとの辞書ファイル。
- `molecular_modeling/structure_out/`, `molecular_modeling/voronoi_pngs/`: 実行後の出力先。
# yoshimi-molecular_modeling
