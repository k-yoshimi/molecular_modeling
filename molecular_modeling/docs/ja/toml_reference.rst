TOML 設定リファレンス
=====================

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
     - グローバル乱数シード。optimize/render の個別シード未設定時のフォールバック
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

置換テーブル
~~~~~~~~~~~~

4 種類の置換テーブルを ``[[テーブル名]]`` 形式で複数定義できます。

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
