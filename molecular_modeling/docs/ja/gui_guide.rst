GUI 操作ガイド
==============

Molecular Modeling GUI は ASE-GUI を拡張したグラフィカルインターフェースです。
分子の追加、置換、エクスポートなどの操作をインタラクティブに行えます。

起動方法
--------

.. code-block:: bash

   # 空の構造で起動
   python molecular_modeling_gui.py

   # 既存の構造ファイルを指定して起動
   python molecular_modeling_gui.py structure.cif

   # CIF ビューアとして起動
   python vis_cif.py structure.cif

   # パイプラインから GUI を起動（edit_and_run.py）
   python edit_and_run.py

ウィンドウ構成
--------------

.. image:: /_static/placeholder.png
   :width: 600
   :alt: GUI ウィンドウ構成

*スクリーンショットをここに差し替えてください*

GUI ウィンドウは以下の要素で構成されます:

1. **メニューバー**: 各種操作メニュー
2. **3D ビューポート**: 分子構造の表示（回転・並進・ズーム可能）
3. **ツールバー**: Translate / Rotate / Zoom ボタン
4. **ステータスバー**: 選択中の原子情報

ウィンドウサイズの変更
~~~~~~~~~~~~~~~~~~~~~~

メニューバーの **Window Size** から変更できます:

- **Small (450x450)**: コンパクト表示
- **Medium (700x700)**: 標準サイズ（デフォルト）
- **Large (1000x800)**: 大画面表示
- **Wide (1200x600)**: 横長表示
- **Custom Size...**: 任意のサイズを指定

ウィンドウの境界をドラッグしてサイズを変更することもできます。

メニュー一覧
------------

MoleculeTools メニュー
~~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - メニュー項目
     - 説明
   * - Load CIF...
     - CIF ファイルを読み込み、構造を表示
   * - Add Molecule...
     - 登録済み分子をシステムに追加
   * - List Molecules/Groups
     - システム中の分子とグループを一覧表示
   * - Save Selection as Molecule...
     - 選択した原子群を新しい分子として TOML に保存
   * - Global Element Substitution...
     - 指定元素を全て別の元素に置換
   * - Molecule Substitution...
     - 特定サイトの分子を別の分子に置換
   * - Molecule Atom Substitution...
     - 分子内の原子リストを新しいリストに置換
   * - Individual Atom Substitution...
     - 1 つの原子を個別に置換

Export メニュー
~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - メニュー項目
     - 説明
   * - Export All (QE, VASP, CIF)...
     - 全フォーマットに一括エクスポート
   * - Quantum ESPRESSO...
     - QE 入力ファイル (.in) としてエクスポート
   * - VASP...
     - POSCAR ファイルとしてエクスポート
   * - CIF...
     - CIF ファイルとしてエクスポート

Window Size メニュー
~~~~~~~~~~~~~~~~~~~~

ウィンドウサイズをプリセットまたはカスタム値で変更できます。

ビュー操作ツールバー
--------------------

3D ビューポートの下部にあるツールバーで、表示を操作できます。
各ボタンは押し続けることで連続的に操作できます。

.. list-table::
   :header-rows: 1
   :widths: 20 30 50

   * - セクション
     - ボタン
     - 説明
   * - Translate
     - ◀ ▶ ▲ ▼
     - ビューを上下左右に移動
   * - Rotate
     - X- X+ Y- Y+ Z- Z+
     - X/Y/Z 軸周りに ±5° 回転
   * - Zoom
     - \- +
     - 拡大・縮小

マウス操作
----------

ASE-GUI 標準のマウス操作が使用できます:

- **左クリック**: 原子を選択
- **Shift + 左クリック**: 選択に追加
- **Ctrl + 左クリック**: 選択を反転
- **左ドラッグ**: ビューの回転
- **右ドラッグ** (または Shift + 右ドラッグ): ビューの並進
- **マウスホイール**: ズーム

ダイアログ詳細
--------------

Add Molecule ダイアログ
~~~~~~~~~~~~~~~~~~~~~~~

.. image:: /_static/placeholder.png
   :width: 400
   :alt: Add Molecule ダイアログ

*スクリーンショットをここに差し替えてください*

- **Molecule**: ドロップダウンから ``data/*.toml`` に定義された分子を選択
- **Position (X, Y, Z)**: 配置位置をデカルト座標で入力
- **Add Mode**:

  - **Single**: 1 つの分子を指定位置に追加
  - **Unit Cell**: 分子定義の全分子を一括追加（格子ベクトルも設定）

- **3D Preview**: 選択した分子のプレビュー表示

Save Selection as Molecule ダイアログ
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. 3D ビューポートで原子を選択（Shift + クリックで複数選択）
2. **MoleculeTools** → **Save Selection as Molecule...** を選択
3. 分子名と単位胞中の分子数を入力
4. **OK** で ``data/<name>.toml`` として保存

Global Element Substitution ダイアログ
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Old Symbol**: 置換元の元素記号（既存元素から選択）
- **New Symbol**: 置換先の元素記号を入力
- 全ての該当原子が一括で置換されます

Molecule Substitution ダイアログ
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Molecule Index**: 置換対象の分子インデックス
  （原子を選択していれば自動設定）
- **New Molecule**: 置換先の分子名をドロップダウンから選択

Export ダイアログ
~~~~~~~~~~~~~~~~~

- 出力先ディレクトリを選択
- 出力フォーマットを選択（QE / VASP / CIF / All）
- **OK** でエクスポート実行
