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
   * - Undo (Ctrl+Z)
     - 直前の操作を取り消す（最大50回）
   * - Redo (Ctrl+Y)
     - 取り消した操作をやり直す
   * - Load Structure... (Ctrl+L)
     - 構造ファイルを読み込み（CIF, XYZ, POSCAR, QE .in, PDB, GEN 等に対応）
   * - Load CIF...
     - CIF ファイルを読み込み、構造を表示
   * - Add Molecule... (Ctrl+J)
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
   * - Batch Element Substitution...
     - 複数の元素置換ルール（例: C→Si, H→F）を一括で定義・実行
   * - Build Supercell...
     - 繰り返し数を指定してスーパーセルを構築
   * - Cell Parameters...
     - 格子定数（a, b, c, α, β, γ）の表示・編集。格子ベクトルと体積も表示。Apply で原子座標をスケーリングして更新
   * - Save Project... (Ctrl+K)
     - プロジェクト（input.toml + data/\*.toml）を指定ディレクトリに保存
   * - Load Project... (Ctrl+T)
     - プロジェクトディレクトリから構造を読み込み
   * - Measurement Log...
     - 測定ログ（距離・角度・二面角の履歴）を表示

Recent Files メニュー
~~~~~~~~~~~~~~~~~~~~

最近開いた CIF ファイルが最大 10 件表示されます。
ファイル名をクリックすると直接開きます。
履歴は ``~/.molecular_modeling_recent.json`` に保存されます。

Export メニュー
~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 35 65

   * - メニュー項目
     - 説明
   * - Export All (QE, VASP, CIF)... (Ctrl+D)
     - 全フォーマットに一括エクスポート
   * - Quantum ESPRESSO...
     - QE 入力ファイル (.in) としてエクスポート
   * - VASP...
     - POSCAR ファイルとしてエクスポート
   * - CIF...
     - CIF ファイルとしてエクスポート

Selection メニュー
~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - メニュー項目
     - 説明
   * - Select by Element...
     - 元素記号を指定して原子を一括選択。複数元素選択可。選択モード（置換/追加）を切り替え可能
   * - Select by Molecule...
     - 分子グループ（mol_id）を指定して原子を一括選択。グループ名・原子数が一覧表示される
   * - Select by Sphere...
     - 選択中の原子を中心に、指定半径内の原子を選択
   * - Delete Selected Atoms
     - 選択した原子を削除（Undo 対応、分子グループメタデータも更新）
   * - Delete Molecule...
     - 選択した原子が属する分子グループ全体を削除
   * - Move Molecule/Atoms...
     - 選択した分子/原子を指定ベクトルで平行移動（分子全体/選択のみ切替可）
   * - Rotate Molecule/Atoms...
     - 選択した分子/原子を重心周りに Euler 角で回転

Pipeline メニュー
~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - メニュー項目
     - 説明
   * - Run Optimize (Voronoi) (Ctrl+P)
     - 現在の構造に対して Voronoi IoU 最適化を実行。結果（IoU 値、Euler 角）がダイアログに表示される
   * - Run Render (Voronoi PNGs)...
     - Voronoi セルオーバーレイの 3D 可視化 PNG を生成
   * - Run Export (QE, VASP, CIF)...
     - 構造ファイルを指定ディレクトリにエクスポート
   * - Run All (Optimize → Render → Export)...
     - 最適化 → 可視化 → エクスポートの全パイプラインを一括実行
   * - Show Optimization Results... (Ctrl+H)
     - 最適化結果（IoU 値、Euler 角、スケール）をテーブル形式で表示。CSV エクスポートおよび Voronoi PNG の閲覧も可能

**注意**: Pipeline メニューの各操作は ``mol_id`` 情報を持つ構造（TOML から構築した構造）に対してのみ実行できます。
処理はバックグラウンドスレッドで実行されるため、GUI は操作中もフリーズしません。

View Settings メニュー
~~~~~~~~~~~~~~~~~~~~~

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - メニュー項目
     - 説明
   * - Atom Appearance...
     - 原子の色・サイズを変更するパネルを表示。色モード（元素別/分子グループ別/タグ別）、元素ごとの色変更、原子サイズスケール、元素別半径の調整が可能
   * - Color by Molecule Group
     - ``mol_id`` に基づき分子グループごとに異なる色で表示
   * - Reset Colors
     - 色設定を ASE デフォルト（jmol）に戻す
   * - Toggle Hydrogen
     - 水素原子の表示/非表示を切り替え
   * - Element Visibility...
     - 元素ごとに表示/非表示を設定するダイアログを表示
   * - Toggle Side Panel
     - サイドパネル（ビュー操作パネル）の表示/非表示を切り替え
   * - Small / Medium / Large / Wide
     - ウィンドウサイズのプリセット変更
   * - Custom Size...
     - 任意のサイズを指定

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
   * - Actions
     - ↩ ↪ ✖ ✥ ↻
     - Undo / Redo / 選択原子の削除 / 移動 / 回転（各ボタンにツールチップ付き）

ツールチップ
~~~~~~~~~~~~

すべてのツールバーボタンにはツールチップ（マウスホバー時のフローティングヘルプ）が表示されます。
ボタンの上にマウスカーソルを置くと、操作内容やキーボードショートカットが黄色のポップアップで表示されます。

サイドパネル
~~~~~~~~~~~~

ウィンドウ右側にサイドパネルが表示されます。よく使うビュー設定にすばやくアクセスできます。

- **Show Hydrogen**: 水素原子の表示/非表示をチェックボックスで切り替え
- **Color Mode**: 色分けモードをラジオボタンで切り替え（元素別 / 分子グループ別 / タグ別）
- **Atom Size**: スライダーで原子サイズスケールを調整（0.3〜3.0）
- **Element Visibility**: 元素ごとに表示/非表示をチェックボックスで個別制御

サイドパネルは **View Settings** → **Toggle Side Panel** で表示/非表示を切り替えられます。

ステータスバー
--------------

画面下部のステータスバーには、選択状態に応じた情報が表示されます。

- **原子未選択時**: 全体情報 — 原子数、含有元素、分子グループ数、セルサイズ
  (例: ``48 atoms | (C, H) | 4 groups | Cell: 5.5x5.4x18.7 Å``)
- **1原子選択時**: 原子詳細 — インデックス、元素名、座標、所属分子名
  (例: ``#3 Carbon (C): 1.114, 0.098, 0.812 Å | mol_id=0 (benzene)``)
- **2原子選択時**: 原子間距離（MIC 補正済み）、所属分子名
  (例: ``#0(C)-#6(H): 1.089 Å | benzene``)
- **3原子選択時**: 結合角度（MIC 補正済み）、所属分子名
  (例: ``#0(C)-#1(C)-#2(C): 120.0° | benzene``)
- **4原子選択時**: 二面角（MIC 補正済み）

測定結果は自動的にログに記録され、 **MoleculeTools** → **Measurement Log...** で一覧を確認できます。

キーボードショートカット
------------------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - ショートカット
     - 操作
   * - Ctrl+Z
     - 元に戻す (Undo)
   * - Ctrl+Y
     - やり直し (Redo)
   * - Ctrl+L
     - 構造ファイルを読み込み（CIF, XYZ, POSCAR 等）
   * - Ctrl+J
     - 分子を追加
   * - Ctrl+D
     - 全フォーマットに一括エクスポート
   * - Ctrl+P
     - Voronoi IoU 最適化を実行
   * - Ctrl+H
     - 最適化結果を表示
   * - Ctrl+K
     - プロジェクトを保存
   * - Ctrl+T
     - プロジェクトを読み込み

ASE-GUI 標準のキーボードショートカット（Ctrl+O で開く、Ctrl+S で保存など）もそのまま使用できます。

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

Atom Appearance パネル
~~~~~~~~~~~~~~~~~~~~~

**View Settings** → **Atom Appearance...** から開くパネルです。3 つのタブがあります:

**Color Mode タブ**:

- **Default (by element)**: 元素ごとの標準色（jmol カラー）
- **By molecule group**: ``mol_id`` に基づき分子グループごとに色分け（20色を循環使用）
- **By tag**: タグ値に応じた色分け
- カスタム配列がある場合、追加の色モードが自動的に表示されます

**Element Colors タブ**:

- 構造中の各元素について個別に色を変更可能
- **Change...** ボタンでカラーピッカーを表示
- **Reset to Default** で jmol デフォルトに戻す

**Size タブ**:

- **Atom Size Scale**: 全体の原子サイズスケール（0.1〜3.0、スライダーで調整）
- **Per-Element Radius Override**: 元素ごとの共有結合半径を個別に変更

Cell Parameters ダイアログ
~~~~~~~~~~~~~~~~~~~~~~~~~

**MoleculeTools** → **Cell Parameters...** から開きます。

- **Lattice lengths**: a, b, c（Å 単位）を表示・編集
- **Lattice angles**: α, β, γ（度単位）を表示・編集
- **Lattice vectors**: 現在の格子ベクトル（読み取り専用）
- **Volume**: 格子体積（Å³）
- **Apply**: 入力値で格子を更新（原子座標はスケーリング）

Export ダイアログ
~~~~~~~~~~~~~~~~~

- 出力先ディレクトリを選択
- 出力フォーマットを選択（QE / VASP / CIF / All）
- **OK** でエクスポート実行
