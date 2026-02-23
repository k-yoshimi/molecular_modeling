入出力ファイル形式
==================

ファイル構成の概要
------------------

本ツールキットは 2 種類の TOML ファイルを使用します:

1. **入力設定** (``input/input.toml``): スーパーセルの格子・繰り返し数・分子配置・置換ルールを定義
2. **分子定義** (``data/<name>.toml``): 各分子の原子構造・格子ベクトル・重心を定義

これらのファイルの関係を以下に示します:

.. code-block:: text

   input/input.toml
   ├── molecules = ["benzene", "naphthalene"]
   │       ↓                    ↓
   │   data/benzene.toml    data/naphthalene.toml
   │   (原子座標・格子)      (原子座標・格子)
   │
   ├── cell = [...]           ← スーパーセルの格子ベクトル
   ├── supercell_repeat = [1,1,2]  ← 繰り返し数
   ├── center_of_mass = [...]      ← 各分子の配置位置
   └── [[substitutions]]           ← 置換ルール（任意）

**処理の流れ**:

1. ``input.toml`` を読み込み、格子ベクトルと繰り返し数を取得
2. ``molecules`` リストに基づき ``data/<name>.toml`` を読み込み
3. ``supercell_repeat`` で指定された繰り返し分のセルを生成
4. 各セルに ``center_of_mass`` の位置で分子を配置
5. 置換ルールを適用（分子置換 → グローバル元素置換 → 分子内原子置換 → 個別原子置換）
6. ASE ``Atoms`` オブジェクトを構築

入力設定: input/input.toml
--------------------------

パイプラインの動作を制御するメインの設定ファイルです。

基本設定
~~~~~~~~

.. code-block:: toml

   # 乱数シード（整数）
   # optimize, render の個別シードが未設定の場合にフォールバック
   rng = 1

   # 格子ベクトル (3×3 配列, 単位: Å)
   # 行で a, b, c ベクトルを定義
   cell = [
     [ 5.522000,  0.000000,  0.000000 ],   # a-ベクトル
     [ 0.000000,  5.439600,  0.000000 ],   # b-ベクトル
     [-3.501254,  0.000000,  9.339672 ],   # c-ベクトル
   ]

   # スーパーセル繰り返し数 [a方向, b方向, c方向]
   supercell_repeat = [1, 1, 2]

   # 使用する分子名のリスト（data/<name>.toml と対応）
   molecules = ["benzene", "naphthalene"]

   # 分率座標を使用するか (true: 分率座標, false: デカルト座標)
   use_fractional = false

   # 各分子の重心位置（molecules と同じ順序）
   center_of_mass = [
     [  0.000000,   0.000000,   0.000000 ],
     [ -1.750627,   2.719800,   4.669836 ],
   ]

全キー一覧
~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - キー
     - 型
     - 説明
   * - ``rng``
     - int
     - グローバル乱数シード。``optimize_rng`` / ``render_rng`` 未設定時のフォールバック
   * - ``optimize_rng``
     - int (任意)
     - Voronoi 最適化ステージ専用の乱数シード。省略時は ``rng`` を使用
   * - ``render_rng``
     - int (任意)
     - 可視化ステージ専用の乱数シード。省略時は ``rng`` を使用
   * - ``data_dir``
     - str (任意)
     - 分子定義 TOML ファイルのディレクトリパス。省略時は ``"./data"``
   * - ``cell``
     - array[3][3]
     - 格子ベクトル（行で a, b, c, 単位 Å）
   * - ``supercell_repeat``
     - array[3]
     - スーパーセル繰り返し数 [nx, ny, nz]
   * - ``molecules``
     - array[str]
     - 使用する分子名（ ``data/<name>.toml`` に対応）
   * - ``use_fractional``
     - bool
     - true の場合 ``center_of_mass`` を分率座標として解釈
   * - ``center_of_mass``
     - array[N][3]
     - 各分子の重心座標（ ``molecules`` と同じ数）

.. _cell-index-explanation-ja:

cell_index とスーパーセルの関係
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``supercell_repeat = [nx, ny, nz]`` でスーパーセルを構築すると、
各セルに ``[ix, iy, iz]`` のインデックスが割り当てられます。

- ``ix``: 0 から nx-1
- ``iy``: 0 から ny-1
- ``iz``: 0 から nz-1

例えば ``supercell_repeat = [1, 1, 2]`` の場合:

.. code-block:: text

   セル [0, 0, 0] — 元の単位胞
   セル [0, 0, 1] — c 軸方向に 1 つずらしたコピー

各セル内に ``molecules`` で定義した分子が配置されるため、
``molecule_index = 0`` はベンゼン、 ``molecule_index = 1`` はナフタレンを指します。

置換テーブルの ``cell_index`` は、このスーパーセル中のセル位置を指定します。

置換テーブル
~~~~~~~~~~~~

4 種類の置換テーブルを ``[[テーブル名]]`` 形式で複数定義できます。
置換は以下の順序で適用されます:

1. **分子置換** (``molecule_substitutions``): 分子全体を別の分子に入れ替え
2. **グローバル元素置換** (``global_atom_substitutions``): 全元素を一括置換
3. **分子内原子置換** (``molecule_atom_substitutions``): 特定分子の原子リストを変更
4. **個別原子置換** (``atom_substitutions``): 1 原子だけ変更

molecule_substitutions（分子置換）
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

特定のセル位置にある分子を別の分子に置き換えます。

.. code-block:: toml

   [[molecule_substitutions]]
   cell_index       = [0, 0, 1]       # スーパーセル中のセル座標
   molecule_index   = 1               # 分子インデックス (0-based)
   new_molecule     = "benzene"       # 置換先の分子名

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - キー
     - 型
     - 説明
   * - ``cell_index``
     - array[3]
     - スーパーセル中のセル位置 [ix, iy, iz]
   * - ``molecule_index``
     - int
     - セル内の分子インデックス (0-based)
   * - ``new_molecule``
     - str
     - 置換先の分子名

global_atom_substitutions（グローバル元素置換）
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

全てのシステムで指定元素を一括置換します。

.. code-block:: toml

   [[global_atom_substitutions]]
   old_symbol = "H"      # 置換元の元素記号
   new_symbol = "He"     # 置換先の元素記号

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - キー
     - 型
     - 説明
   * - ``old_symbol``
     - str
     - 置換元の元素記号
   * - ``new_symbol``
     - str
     - 置換先の元素記号

molecule_atom_substitutions（分子内原子置換）
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

特定の分子インスタンスの原子リストを新しいリストに置き換えます。

.. note::

   この置換は原子リスト全体を上書きするため、同じ分子に対する
   グローバル元素置換（``global_atom_substitutions``）の結果は無効になります。

.. code-block:: toml

   [[molecule_atom_substitutions]]
   cell_index     = [0, 0, 1]
   molecule_index = 1
   new_symbols    = ["Si6O6"]    # コンパクト表記可

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - キー
     - 型
     - 説明
   * - ``cell_index``
     - array[3]
     - スーパーセル中のセル位置
   * - ``molecule_index``
     - int
     - セル内の分子インデックス
   * - ``new_symbols``
     - array[str]
     - 新しい元素記号リスト（コンパクト表記可）

atom_substitutions（個別原子置換）
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

特定の 1 原子だけを置き換えます。

.. code-block:: toml

   [[atom_substitutions]]
   cell_index     = [0, 0, 0]
   molecule_index = 1
   atom_index     = 2           # 分子内の原子インデックス (0-based)
   new_symbol     = "Fe"

.. list-table::
   :header-rows: 1
   :widths: 25 15 60

   * - キー
     - 型
     - 説明
   * - ``cell_index``
     - array[3]
     - スーパーセル中のセル位置
   * - ``molecule_index``
     - int
     - セル内の分子インデックス
   * - ``atom_index``
     - int
     - 分子内の原子インデックス (0-based)
   * - ``new_symbol``
     - str
     - 置換先の元素記号

分子定義: data/<name>.toml
--------------------------

各分子の構造情報を定義するファイルです。
``input.toml`` の ``molecules`` リストに記載された名前と
ファイル名（ ``.toml`` 拡張子を除く）が対応します。

.. note::

   ``data/*.toml`` の ``center_of_mass`` は **分子辞書の参照情報** として使用されます。
   スーパーセル中での分子配置位置は ``input.toml`` の ``center_of_mass`` で決定されます。

   ``data/*.toml`` の ``lattice_vectors`` は **分子固有の格子** です。
   スーパーセル全体の格子は ``input.toml`` の ``cell`` で定義されます。
   混合結晶（異なる分子の共結晶）の場合、これらは異なる値になります。

基本構造
~~~~~~~~

.. code-block:: toml

   molecule_name = "benzene"

   # 元素記号: 個別リストまたはコンパクト表記
   symbols = ["C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H"]
   # または: symbols = ["C6H6"]

   # 座標系の指定
   use_fractional = false

   # 原子座標 (N×3 配列)
   positions = [
     [-1.113710, -0.097913, -0.811833],
     [ 1.113710,  0.097913,  0.811833],
     # ...
   ]

   # 格子ベクトル (3×3 配列, 単位: Å)
   lattice_vectors = [
     [ 5.522000,  0.000000,  0.000000],
     [ 0.000000,  5.439600,  0.000000],
     [-2.693272,  0.000000,  7.184363],
   ]

   # 単位胞中の分子数
   num_molecules = 2

   # 各分子の重心位置
   center_of_mass = [
     [ 0.000000,  0.000000,  0.000000],
     [-1.346636,  2.719800,  3.592182],
   ]

   # Voronoi 参照分子インデックス (省略時: 0)
   # voronoi_ref_id = 0

全キー一覧
~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 25 10 10 55

   * - キー
     - 型
     - 必須
     - 説明
   * - ``molecule_name``
     - str
     - Yes
     - 分子の識別名
   * - ``symbols``
     - array[str]
     - Yes
     - 元素記号のリスト。コンパクト表記可（例: ``"C6H6"``）
   * - ``use_fractional``
     - bool
     - No
     - true の場合、positions と center_of_mass を分率座標として解釈
   * - ``positions``
     - array[N][3]
     - Yes
     - 原子の座標（デカルト Å または分率座標）
   * - ``lattice_vectors``
     - array[3][3]
     - Yes
     - 格子ベクトル（行で a, b, c）
   * - ``num_molecules``
     - int
     - Yes
     - 単位胞中の分子数
   * - ``center_of_mass``
     - array[K][3]
     - Yes
     - K 個の分子の重心座標
   * - ``voronoi_ref_id``
     - int
     - No
     - Voronoi セル計算の参照分子インデックス (デフォルト: 0)

コンパクト表記
~~~~~~~~~~~~~~

``symbols`` フィールドでは、元素記号と数を連結した表記が使えます:

.. code-block:: text

   "C6H6"    → ["C", "C", "C", "C", "C", "C", "H", "H", "H", "H", "H", "H"]
   "C10H8"   → ["C"] × 10 + ["H"] × 8
   "Si6O6"   → ["Si"] × 6 + ["O"] × 6

座標系の変換
~~~~~~~~~~~~

``use_fractional = true`` の場合、内部で以下の変換が行われます:

.. code-block:: text

   pos_cart = pos_frac @ lattice_vectors
   pos_frac = pos_cart @ inv(lattice_vectors)

.. note::

   ``input.toml`` と ``data/*.toml`` の ``use_fractional`` は **独立** です。
   一方が分率座標を使い、もう一方がデカルト座標を使うことも可能です。

完全な設定例
------------

以下に、同梱のベンゼン・ナフタレン混合結晶の完全な設定例を示します。

input/input.toml の完全例
~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: toml

   # 乱数シード
   rng = 1

   # スーパーセルの格子ベクトル (Å)
   # ベンゼンとナフタレンの混合結晶の格子
   cell = [
     [ 5.522000,  0.000000,  0.000000 ],
     [ 0.000000,  5.439600,  0.000000 ],
     [-3.501254,  0.000000,  9.339672 ],
   ]

   # c 軸方向に 2 倍のスーパーセルを構築
   supercell_repeat = [1, 1, 2]

   # 2 種類の分子を使用
   # data/benzene.toml と data/naphthalene.toml を参照
   molecules = ["benzene", "naphthalene"]

   # デカルト座標を使用
   use_fractional = false

   # 各分子の配置位置（molecules と同じ順序）
   # benzene → 原点, naphthalene → オフセット位置
   center_of_mass = [
     [  0.000000,   0.000000,   0.000000 ],
     [ -1.750627,   2.719800,   4.669836 ],
   ]

   # --- 以下はオプションの置換ルール ---

   # セル [0,0,1] の 2 番目の分子(naphthalene)を benzene に置換
   [[molecule_substitutions]]
   cell_index       = [0, 0, 1]
   molecule_index   = 1
   new_molecule     = "benzene"

   # 全ての H を He に置換
   [[global_atom_substitutions]]
   old_symbol = "H"
   new_symbol = "He"

   # セル [0,0,1] の 2 番目の分子の原子リストを Si6O6 に変更
   [[molecule_atom_substitutions]]
   cell_index     = [0, 0, 1]
   molecule_index = 1
   new_symbols    = ["Si6O6"]

   # セル [0,0,0] の 2 番目の分子の 3 番目の原子を Fe に変更
   [[atom_substitutions]]
   cell_index     = [0, 0, 0]
   molecule_index = 1
   atom_index     = 2
   new_symbol     = "Fe"

data/benzene.toml の完全例
~~~~~~~~~~~~~~~~~~~~~~~~~~

元素記号を **個別リスト** で記述する形式です。

.. code-block:: toml

   molecule_name = "benzene"

   # 12 原子を個別に列挙
   symbols = [
     "C", "C", "C", "C", "C", "C",
     "H", "H", "H", "H", "H", "H",
   ]

   use_fractional = false

   # 第 1 分子の原子座標 (Å)
   # symbols と同じ順序で対応
   positions = [
     # 炭素原子 (6 個)
     [ -1.113710, -0.097913, -0.811833 ],
     [  1.113710,  0.097913,  0.811833 ],
     [ -0.115564, -1.087920, -0.865716 ],
     [  0.115564,  1.087920,  0.865716 ],
     [ -0.979384,  0.951930,  0.052446 ],
     [  0.979384, -0.951930, -0.052446 ],
     # 水素原子 (6 個)
     [ -1.861055, -0.158292, -1.362155 ],
     [  1.861055,  0.158292,  1.362155 ],
     [ -0.197159, -1.816826, -1.437591 ],
     [  0.197159,  1.816826,  1.437591 ],
     [ -1.651612,  1.593259,  0.089805 ],
     [  1.651612, -1.593259, -0.089805 ],
   ]

   # ベンゼン結晶の格子ベクトル
   lattice_vectors = [
     [ 5.522000,  0.000000,  0.000000],
     [ 0.000000,  5.439600,  0.000000],
     [-2.693272,  0.000000,  7.184363],
   ]

   # 単位胞中に 2 分子
   num_molecules = 2

   # 各分子の重心位置
   center_of_mass = [
     [ 0.000000,  0.000000,  0.000000],
     [-1.346636,  2.719800,  3.592182],
   ]

data/naphthalene.toml の完全例
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

元素記号を **コンパクト表記** で記述する形式です。

.. code-block:: toml

   molecule_name = "naphthalene"

   # コンパクト表記: C10H8 → C×10 + H×8 = 18 原子
   symbols = ["C10H8"]

   use_fractional = false

   # 第 1 分子の原子座標 (Å)
   # コンパクト表記の展開順 (C が先、H が後) と対応
   positions = [
     # 炭素原子 (10 個)
     [ -0.873883,  0.110730,  2.365899 ],
     [  0.873883, -0.110730, -2.365899 ],
     [ -0.125794,  0.966533,  1.590346 ],
     [  0.125794, -0.966533, -1.590346 ],
     [  0.213579,  0.617262,  0.251858 ],
     [ -0.213579, -0.617262, -0.251858 ],
     [  0.981797,  1.483077, -0.575778 ],
     [ -0.981797, -1.483077,  0.575778 ],
     [ -1.309714, -1.131450,  1.848626 ],
     [  1.309714,  1.131450, -1.848626 ],
     # 水素原子 (8 個)
     [ -1.094109,  0.351038,  3.237058 ],
     [  1.094109, -0.351038, -3.237058 ],
     [  0.159973,  1.778750,  1.942805 ],
     [ -0.159973, -1.778750, -1.942805 ],
     [  1.263671,  2.303540, -0.240443 ],
     [ -1.263671, -2.303540,  0.240443 ],
     [ -1.818179, -1.705126,  2.375888 ],
     [  1.818179,  1.705126, -2.375888 ],
   ]

   # ナフタレン結晶の格子ベクトル（ベンゼンとは異なる）
   lattice_vectors = [
     [ 8.034800,  0.000000,  0.000000],
     [ 0.000000,  5.889900,  0.000000],
     [-4.738553,  0.000000,  7.134798],
   ]

   # 単位胞中に 2 分子
   num_molecules = 2

   # 各分子の重心位置
   center_of_mass = [
     [ 0.000000,  0.000000,  0.000000],
     [-1.346636,  2.719800,  3.592182],
   ]

新しい分子定義の作成方法
~~~~~~~~~~~~~~~~~~~~~~~~

独自の分子を追加するには、以下の手順で ``data/<name>.toml`` を作成します:

1. ``molecule_name`` に分子名を設定
2. ``symbols`` で元素記号を定義（個別リストまたはコンパクト表記）
3. ``positions`` に **第 1 分子** の原子座標を記述（重心が原点付近になるように）
4. ``lattice_vectors`` にその分子の結晶格子を設定
5. ``num_molecules`` に単位胞あたりの分子数を設定
6. ``center_of_mass`` に各分子の重心位置を記述

.. tip::

   GUI の **MoleculeTools** → **Save Selection as Molecule...** を使うと、
   3D ビューで選択した原子群を ``data/*.toml`` として保存できます。
