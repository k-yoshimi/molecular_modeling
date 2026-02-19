チュートリアル
==============

このチュートリアルでは、ベンゼンとナフタレンの混合結晶を例に、
パイプラインの全ステップを実行します。

ステップ 1: 分子定義の理解
--------------------------

まず、 ``data/benzene.toml`` の内容を見てみましょう:

.. code-block:: toml

   molecule_name = "benzene"

   symbols = [
     "C", "C", "C", "C", "C", "C",
     "H", "H", "H", "H", "H", "H",
   ]

   use_fractional = false

   positions = [
     [ -1.113710, -0.097913, -0.811833 ],
     [  1.113710,  0.097913,  0.811833 ],
     # ... (12原子分の座標)
   ]

   lattice_vectors = [
     [ 5.522000,  0.000000,  0.000000],
     [ 0.000000,  5.439600,  0.000000],
     [-2.693272,  0.000000,  7.184363],
   ]

   num_molecules = 2

   center_of_mass = [
     [ 0.000000,  0.000000,  0.000000],
     [-1.346636,  2.719800,  3.592182],
   ]

各フィールドの意味:

- ``molecule_name``: 分子名（ ``input.toml`` の ``molecules`` で参照）
- ``symbols``: 元素記号のリスト。 ``"C6H6"`` のようなコンパクト表記も可
- ``positions``: 原子のデカルト座標 (Å)
- ``lattice_vectors``: 格子ベクトル（行で a, b, c を定義）
- ``num_molecules``: 単位胞中の分子数
- ``center_of_mass``: 各分子の重心位置

ナフタレン（ ``data/naphthalene.toml`` ）も同様の構造です。
元素記号はコンパクト表記 ``"C10H8"`` を使用しています。

ステップ 2: 入力ファイルの設定
------------------------------

``input/input.toml`` でスーパーセルの構築方法を指定します:

.. code-block:: toml

   # 乱数シード
   rng = 1

   # 格子ベクトル (Å)
   cell = [
     [ 5.522000,  0.000000,  0.000000 ],
     [ 0.000000,  5.439600,  0.000000 ],
     [-3.501254,  0.000000,  9.339672 ],
   ]

   # スーパーセルの繰り返し数 [a, b, c方向]
   supercell_repeat = [1, 1, 2]

   # 使用する分子名のリスト
   molecules = ["benzene", "naphthalene"]

   # 分率座標を使用するか
   use_fractional = false

   # 各分子の重心位置
   center_of_mass = [
     [  0.000000,   0.000000,   0.000000 ],
     [ -1.750627,   2.719800,   4.669836 ],
   ]

``supercell_repeat = [1, 1, 2]`` は c 軸方向に 2 倍のスーパーセルを作ることを意味します。
2 種類の分子（ベンゼンとナフタレン）がそれぞれの重心位置に配置されます。

ステップ 3: パイプラインの実行
------------------------------

.. code-block:: bash

   python run.py

パイプラインは以下の順序で処理を行います:

1. **Build**: TOML を読み込み、ASE ``Atoms`` オブジェクトを構築
2. **Group**: 原子を ``mol_id`` でグループ化
3. **Optimize**: Euler 角のグリッドサーチで Voronoi IoU を最大化
4. **Check**: Voronoi 体積のタイリング整合性を検証
5. **Render**: Voronoi セルオーバーレイの 3D 可視化
6. **Export**: QE, VASP, CIF ファイルを出力

ステップ 4: 出力ファイルの確認
------------------------------

構造ファイル
~~~~~~~~~~~~

``source/structure_out/`` に以下のファイルが生成されます:

**qe.in** (Quantum ESPRESSO):

.. code-block:: text

   &CONTROL
     calculation = 'scf'
     prefix = 'molecule'
     outdir = './tmp'
   /
   &SYSTEM
     ecutwfc = 40
     ecutrho = 320
   /
   ...

**POSCAR** (VASP): 原子位置が分率座標で記述されます。

**structure.cif**: 汎用 CIF フォーマットで出力されます。

ステップ 5: Voronoi 可視化の見方
--------------------------------

``source/voronoi_pngs/`` に生成される PNG 画像は、各サイトの Voronoi セル
オーバーレイを示しています:

- **青い半透明ポリヘドロン**: システムの Voronoi セル（実際のセル分割）
- **赤い半透明ポリヘドロン**: 辞書の Voronoi セル（最適化ターゲット）
- IoU (Intersection over Union) 値が高いほど、分子配向の再現性が高い

ステップ 6: GUI でのインタラクティブ操作
----------------------------------------

CIF ファイルの読み込み
~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   python vis_cif.py source/structure_out/structure.cif

または GUI 内のメニューから:

**MoleculeTools** → **Load CIF...**

.. image:: /_static/placeholder.png
   :width: 400
   :alt: GUI メインウィンドウ

*スクリーンショットをここに差し替えてください*

分子の追加
~~~~~~~~~~

1. **MoleculeTools** → **Add Molecule...** を選択
2. ドロップダウンから分子名を選択（例: benzene）
3. 配置位置の座標を入力
4. **OK** をクリック

置換操作
~~~~~~~~

GUI では 4 種類の置換操作が可能です:

1. **Global Element Substitution**: 全ての指定元素を一括置換
2. **Molecule Substitution**: 特定サイトの分子を別の分子に置換
3. **Molecule Atom Substitution**: 分子内の原子リストを変更
4. **Individual Atom Substitution**: 個別の原子を置換

エクスポート
~~~~~~~~~~~~

**Export** → **Export All (QE, VASP, CIF)...** で全フォーマットに一括エクスポートできます。

ステップ 7: 置換機能の詳細
--------------------------

``input.toml`` では 4 種類の置換テーブルを設定できます。

分子置換 (molecule_substitutions)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

特定のセル位置にある分子を別の分子に置き換えます:

.. code-block:: toml

   [[molecule_substitutions]]
   cell_index       = [0, 0, 1]    # スーパーセル中のセル位置
   molecule_index   = 1            # 分子インデックス (0-based)
   new_molecule     = "benzene"    # 置換先の分子名

グローバル元素置換 (global_atom_substitutions)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

全ての指定元素を一括で置き換えます:

.. code-block:: toml

   [[global_atom_substitutions]]
   old_symbol = "H"     # 置換元の元素
   new_symbol = "He"    # 置換先の元素

分子内原子置換 (molecule_atom_substitutions)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

特定の分子インスタンスの全原子を新しいリストで置き換えます:

.. code-block:: toml

   [[molecule_atom_substitutions]]
   cell_index     = [0, 0, 1]
   molecule_index = 1
   new_symbols    = ["Si6O6"]    # コンパクト表記可

個別原子置換 (atom_substitutions)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1 つの原子だけを置き換えます:

.. code-block:: toml

   [[atom_substitutions]]
   cell_index     = [0, 0, 0]
   molecule_index = 1
   atom_index     = 2           # 分子内の原子インデックス (0-based)
   new_symbol     = "Fe"
