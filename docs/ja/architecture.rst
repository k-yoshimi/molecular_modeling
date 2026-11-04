アーキテクチャ
==============

パイプラインフロー
------------------

メインパイプライン（ ``run.py`` → ``source/main.py`` ）は以下の順序で処理を実行します:

.. code-block:: text

   input/input.toml + data/*.toml
          │
          ▼
   ┌─────────────────────────┐
   │  build_system_from_toml │  TOML → ASE Atoms 構築
   │                         │  スーパーセル生成、置換適用
   └──────────┬──────────────┘
              │  system (Atoms + system.info)
              ▼
   ┌─────────────────────────┐
   │  build_molecule_groups  │  mol_id による原子グループ化
   └──────────┬──────────────┘
              │  groups: list[list[int]]
              ▼
   ┌─────────────────────────┐
   │ optimize_voronoi_overlap│  Euler 角グリッドサーチ +
   │                         │  Nelder-Mead 最適化
   └──────────┬──────────────┘
              │  optimized_angles_deg, scales, iou
              ▼
   ┌──────────────────────────────────┐
   │ check_system_voronoi_volume     │  Voronoi 体積の
   │ _tiling (voronoi_geom.py)       │  タイリング整合性検証
   └──────────┬───────────────────────┘
              │
              ▼
   ┌─────────────────────────┐
   │ render_voronoi_overlays │  3D matplotlib 可視化
   └──────────┬──────────────┘
              │  voronoi_pngs/*.png
              ▼
   ┌─────────────────────────┐
   │     write_structure     │  QE / VASP / CIF エクスポート
   └──────────┬──────────────┘
              │
              ▼
   structure_out/  (qe.in, POSCAR, INCAR, KPOINTS, structure.cif)

各モジュールの役割
------------------

build_system_from_toml.py
~~~~~~~~~~~~~~~~~~~~~~~~~

**入力**: ``input/input.toml``, ``data/*.toml``

**出力**: ASE ``Atoms`` オブジェクト

主な処理:

1. TOML ファイルの読み込みとパース
2. 分率座標 → デカルト座標への変換（ ``use_fractional`` が true の場合）
3. スーパーセルの構築（繰り返しタイリング）
4. 各分子の配置（重心位置への並進）
5. 4 種類の置換テーブルの適用
6. ``mol_id`` 配列の生成（各原子がどのサイトに属するか）
7. ``system.info`` へのメタデータ設定

build_molecule_groups.py
~~~~~~~~~~~~~~~~~~~~~~~~

**入力**: ``Atoms`` オブジェクト（ ``mol_id`` 配列付き）

**出力**: ``list[list[int]]`` — サイトごとの原子インデックスリスト

``mol_id`` 配列を走査して、同じ値を持つ原子のインデックスをグループ化します。

voronoi_geom.py
~~~~~~~~~~~~~~~

コアのジオメトリモジュールです。以下を提供します:

- **Polyhedron クラス**: 凸多面体の表現（頂点 + 面）
- **Voronoi セル計算**: SciPy の halfspace intersection による周期的 Voronoi 分割
- **IoU 計算**: trimesh と manifold3d を使った体積ベースの Intersection over Union
- **回転・スケーリング**: Euler 角による多面体の回転と均一スケーリング
- **周期境界処理**: セル境界でのラッピング

optimize_voronoi_overlap.py
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**処理**: サイトごとに Voronoi IoU を最大化する Euler 角を探索

アルゴリズム:

1. システム Voronoi セル体積 ``V_sys`` を計算
2. 辞書 Voronoi セル体積 ``V_dict`` を計算
3. 均一スケール ``s = (V_sys / V_dict)^(1/3)`` を設定
4. Euler 角 (α, β, γ) の粗いグリッドサーチ（デフォルト 13×13×13）
5. 最良点から Nelder-Mead 法で精密最適化
6. 最適な回転とスケーリングを原子位置に適用

結果は ``system.info`` に格納:

- ``optimized_angles_deg``: 最適 Euler 角 (N×3)
- ``optimized_scales``: 最適スケール (N,)
- ``optimized_iou``: IoU 値 (N,)

render_voronoi_overlays.py
~~~~~~~~~~~~~~~~~~~~~~~~~~

matplotlib の 3D プロットでシステムと辞書の Voronoi セルを重ね合わせて可視化します。

- 青い半透明: システム Voronoi セル
- 赤い半透明: 辞書 Voronoi セル
- サイトごとに 1 枚の PNG を生成

write_structure.py
~~~~~~~~~~~~~~~~~~

ASE の I/O 機能を使って、以下のフォーマットで構造ファイルを出力します:

- **Quantum ESPRESSO**: ``qe.in`` (SCF 計算設定付き)
- **VASP**: ``POSCAR``, ``INCAR``, ``KPOINTS``
- **CIF**: ``structure.cif``

molecule_unwrap.py
~~~~~~~~~~~~~~~~~~

周期境界で切れた分子をつなぎ直します（結合グラフの幅優先探索。結合の判定は自然共有結合半径の 1.2 倍）。
原子の集まりがつながった 1 分子かどうかの判定と、結合している原子対の列挙も行います。
GUI（Save Selection as Molecule）と CSP パッケージから使います。

csp/ — 結晶構造予測
~~~~~~~~~~~~~~~~~~~

孤立分子 1 個から結晶の候補を作る別のパイプラインです
（``molecular-modeling-csp`` → ``source/csp/cli.py`` → ``pipeline.run_csp``）。
``csp`` extra（PyXtal、tblite、pymatgen）が必要です。:doc:`csp_guide` を参照してください。

.. code-block:: text

   input/csp.toml + 分子ファイル
          |
          v
   候補生成 (gen_random, gen_voronoi, gen_seed) -> 事前ふるい (screen)
          -> 緩和 (relax, GFN2-xTB) -> 重複除去 (dedupe)
          -> csp_out/ (candidates.csv, NNN_<id>.cif, run.log)

- ``config.py``: ``CSPConfig`` と ``load_csp_config``\ （TOML 設定、blind/informed モード、値の範囲チェック）
- ``pipeline.py``: ``run_csp``\ （各段、時間の目安、依存ライブラリの確認、出力）。``cli.py``: コマンドラインの入口
- ``molecule.py``: 分子の読み込み、体積、RMSD、剛体の重ね合わせ（``match_and_align``）
- ``gen_random.py``・``gen_voronoi.py``・``gen_seed.py``: ランダム（PyXtal）、Voronoi による向き決め、置換元からの生成器
- ``screen.py``: 重なりと充填率のチェック、剛体 Lennard-Jones による事前ふるい
- ``relax.py``: 並列の緩和と、変形・結合のチェック
- ``dedupe.py``: 同じ構造をまとめる
- ``candidate.py``: 段の間で受け渡す ``Candidate``
- ``deps.py``: ``csp`` extra が import できるかの確認
- ``outdir.py``: 実行と答え合わせが共有するロックファイルと、結果が未完成であることを示す目印
- ``validate.py``: 実験 CIF との比較（``python -m molecular_modeling.csp.validate``）

system.info メタデータ
----------------------

パイプラインの各ステージ間でデータを受け渡すために、
ASE ``Atoms`` オブジェクトの ``info`` 辞書を使用します。

.. list-table::
   :header-rows: 1
   :widths: 30 15 55

   * - キー
     - 型
     - 説明
   * - ``groups``
     - list[list[int]]
     - サイトごとの原子インデックス
   * - ``site_mol_name``
     - list[str]
     - サイトごとの分子名
   * - ``dict_meta``
     - dict
     - 分子ごとのメタデータ（格子ベクトル、重心、voronoi_ref_id）
   * - ``rng``
     - dict
     - 乱数シード（global, optimize, render）
   * - ``optimized_angles_deg``
     - ndarray(N,3)
     - 最適化後の Euler 角（度）
   * - ``optimized_scales``
     - ndarray(N,)
     - 最適化後のスケール係数
   * - ``optimized_iou``
     - ndarray(N,)
     - 最適化後の IoU 値

原子配列
~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 20 15 65

   * - 配列名
     - 型
     - 説明
   * - ``mol_id``
     - int[N_atoms]
     - 各原子が属するサイトのインデックス

GUI アーキテクチャ
------------------

GUI は MVC (Model-View-Controller) パターンに近い設計です。

.. code-block:: text

   ┌────────────────────────────────────────────────────────────┐
   │              MolecularModelingGUI (View)                  │
   │                                                            │
   │  メニューバー: MoleculeTools / Recent Files / Export /     │
   │              Selection / Pipeline / View Settings / Help   │
   │                                                            │
   │  ┌─────────────────┐  ┌──────────────────────────────┐   │
   │  │  ファイルパネル  │  │  サイドパネル (右)            │   │
   │  │  (Input/Output/ │  │  - カラーモード               │   │
   │  │   Edit/Pipeline)│  │  - ボンド表示 / 3D 球体      │   │
   │  └─────────────────┘  │  - 元素コントロール           │   │
   │  ┌─────────────────┐  │  - 分子コントロール           │   │
   │  │  3D Viewport    │  │  - 原子リスト                 │   │
   │  │  (ASE Canvas)   │  └──────────────────────────────┘   │
   │  └─────────────────┘                                      │
   │  ┌─────────────────┐  ┌──────────────────────────────┐   │
   │  │  ツールバー      │  │  ボトムパネル                 │   │
   │  │  (移動/回転/     │  │  (Molecule Tools /            │   │
   │  │   ズーム/BG色)   │  │   Structure Tools)            │   │
   │  └─────────────────┘  └──────────────────────────────┘   │
   │  ステータスバー (選択情報・距離・角度・二面角)              │
   └────────────────────────────────────┬───────────────────────┘
                                        │ ダイアログ呼び出し
                                        ▼
   ┌────────────────────────────────────────────────────────────┐
   │                     Dialogs (View)                        │
   │  AddMolecule / Substitution (4種) / BatchSubstitution     │
   │  Export / SaveMolecule / MoleculeList / Supercell         │
   │  CellParameters / Focus / Move / Rotate / SelectBy...     │
   └────────────────────────────┬───────────────────────────────┘
                                │ 結果を返す
                                ▼
   ┌────────────────────────────────────────────────────────────┐
   │           MolecularController (Controller)                │
   │  add_molecule / substitute_* / build_supercell /          │
   │  run_optimize / run_render / run_export / get_atoms       │
   └────────────────────────────┬───────────────────────────────┘
                                │ Atoms を操作
                                ▼
   ┌────────────────────────────────────────────────────────────┐
   │               ASE Atoms (Model)                           │
   │  positions / symbols / cell / pbc / info / arrays         │
   └────────────────────────────────────────────────────────────┘

GUI モジュール構成:

- ``molecular_modeling_gui.py``: メインの GUI クラス（ASE-GUI を拡張）。ウィンドウレイアウト、メニュー、ツールバー、サイドパネル、ステータスバーを管理
- ``source/gui/controller.py``: GUI 操作をバックエンドに橋渡し
- ``source/gui/theme.py``: テーマ設定（配色・フォント・余白の定数、ttk テーマ構成）
- ``source/gui/base_dialog.py``: モーダルダイアログの基底クラス（OK/Cancel, バリデーション, 座標入力）
- ``source/gui/molecule_registry.py``: ``data/*.toml`` から分子定義を動的に読み込み
- ``source/gui/add_molecule_dialog.py``: 分子追加ダイアログ
- ``source/gui/substitution_dialogs.py``: 4 種類の置換ダイアログ + 一括元素置換ダイアログ
- ``source/gui/export_dialog.py``: エクスポートダイアログ
- ``source/gui/save_molecule_dialog.py``: 選択原子の分子保存ダイアログ
- ``source/gui/molecule_list_dialog.py``: 分子一覧ダイアログ
- ``source/gui/supercell_dialog.py``: スーパーセル設定ダイアログ
