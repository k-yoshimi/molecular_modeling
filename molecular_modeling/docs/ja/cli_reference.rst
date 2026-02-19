CLI リファレンス
================

run.py — パイプライン実行
-------------------------

メインのパイプラインを実行します。

.. code-block:: bash

   python run.py

**処理フロー**:

1. ``input/input.toml`` と ``data/*.toml`` からシステムを構築
2. 原子を分子グループにインデックス化
3. Voronoi IoU 最大化による分子配向の最適化
4. Voronoi 体積のタイリング整合性チェック
5. Voronoi セルオーバーレイの可視化 → ``source/voronoi_pngs/``
6. 構造ファイルの出力 → ``source/structure_out/``

**入力ファイル**:

- ``input/input.toml`` — 主設定
- ``data/<name>.toml`` — 分子定義

**出力ファイル**:

- ``source/structure_out/qe.in`` — Quantum ESPRESSO 入力
- ``source/structure_out/POSCAR`` — VASP 構造ファイル
- ``source/structure_out/INCAR`` — VASP 設定ファイル
- ``source/structure_out/KPOINTS`` — VASP k点設定
- ``source/structure_out/structure.cif`` — CIF ファイル
- ``source/voronoi_pngs/*.png`` — Voronoi 可視化画像

molecular_modeling_gui.py — GUI 起動
------------------------------------

ASE-GUI を拡張した分子モデリング GUI を起動します。

.. code-block:: bash

   # 空の構造で起動
   python molecular_modeling_gui.py

   # 構造ファイルを指定して起動
   python molecular_modeling_gui.py [structure_file]

**引数**:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - 引数
     - 説明
   * - ``structure_file`` (任意)
     - ASE が読み込み可能な構造ファイルのパス

ファイルを指定しない場合、20×20×20 Å のボックスが生成されます。

vis_cif.py — CIF ビューア
--------------------------

CIF ファイルを GUI で表示します。

.. code-block:: bash

   # ファイルを指定して起動
   python vis_cif.py <path_to_file.cif>

   # ファイルダイアログで選択
   python vis_cif.py

**引数**:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - 引数
     - 説明
   * - ``path_to_file.cif`` (任意)
     - CIF ファイルのパス。省略時はファイル選択ダイアログが表示される

**備考**: ASE の CIF パーサーで読み込みに失敗した場合、
``_citation*`` タグを除去して再試行するフォールバック機能があります。

edit_and_run.py — GUI 編集＋最適化
-----------------------------------

GUI での手動編集後に Voronoi 最適化を実行するパイプラインです。

.. code-block:: bash

   python edit_and_run.py

**処理フロー**:

1. TOML からシステムを構築
2. **GUI を起動** （ユーザーが構造を編集）
3. GUI を閉じると、編集後のシステムで処理を続行
4. 分子グループの構築
5. Voronoi IoU 最適化
6. 可視化とエクスポート

テストの実行
------------

.. code-block:: bash

   # 全テスト実行
   pytest

   # 特定のテストファイル
   pytest tests/test_voronoi_geom.py

   # 特定のテストケース
   pytest tests/test_voronoi_geom.py::TestPolyhedron::test_volume

   # 遅いテストをスキップ
   pytest -m "not slow"
