# molecular_modeling

Python/ASE による分子配置最適化と構造ファイル生成スクリプトの実体です。`input/input.toml` と `data/*.toml` をもとにスーパーセルを組み立て、Voronoi 重なり最小化 → 画像出力 → 各種計算コード向け入力生成までを行います。

## 使い方
- 実行場所: このディレクトリ直下。
- コマンド:
  ```bash
  python run.py
  ```
- 主な処理:
  1. `source/build_system_from_toml.py` で入力 TOML から ASE `Atoms` を構築。
  2. `source/optimize_voronoi_overlap.py` で分子配向を最適化。
  3. `source/render_voronoi_overlays.py` で Voronoi セル重ね合わせ画像を `voronoi_pngs/` へ保存。
  4. `source/write_structure.py` で `structure_out/` に QE/ VASP/ CIF を書き出し。

## 入力設定の要点
- `input/input.toml`: セルベクトル、`supercell_repeat`、`molecules`、`center_of_mass` などを指定。`use_fractional=true` で分数座標解釈に切替可。
- 置換: `molecule_substitutions`, `atom_substitutions`, `global_atom_substitutions`, `molecule_atom_substitutions` で分子/原子種を差し替え。
- 乱数: `rng`（または `optimize_rng` / `render_rng`）でシードを制御。
- 分子辞書: `data/<name>.toml` に `positions`, `center_of_mass`, `lattice_vectors`, 追加で `voronoi_ref_id` などを定義。

## 出力
- `structure_out/`: `qe.in`, `POSCAR`, `INCAR`, `KPOINTS`, `structure.cif`, `POSCAR.vasp` など。
- `voronoi_pngs/`: 最適化後の Voronoi 可視化 PNG。

## 主な依存モジュール
- コア計算: `ase`, `numpy`, `scipy`
- 可視化: `matplotlib`, `contourpy`, `pillow`
- 形状処理: `trimesh`, `manifold3d`
- その他: `packaging`, `python-dateutil` など（詳細は `../requirements.txt`）
