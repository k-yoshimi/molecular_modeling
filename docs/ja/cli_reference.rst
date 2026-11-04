コマンドラインツール
====================

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
5. Voronoi セルオーバーレイの可視化 → ``voronoi_pngs/``
6. 構造ファイルの出力 → ``structure_out/``

**入力ファイル** （詳細は :doc:`toml_reference` を参照）:

- ``input/input.toml`` — 主設定
- ``data/<name>.toml`` — 分子定義

**出力ファイル**:

- ``structure_out/qe.in`` — Quantum ESPRESSO 入力
- ``structure_out/POSCAR`` — VASP 構造ファイル
- ``structure_out/INCAR`` — VASP 設定ファイル
- ``structure_out/KPOINTS`` — VASP k点設定
- ``structure_out/structure.cif`` — CIF ファイル
- ``voronoi_pngs/*.png`` — Voronoi 可視化画像

molecular_modeling_gui.py — GUI 起動
------------------------------------

ASE-GUI を拡張した分子モデリング GUI を起動します。
詳細は :doc:`gui_guide` を参照してください。

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

molecular-modeling-csp — 結晶構造予測
-------------------------------------

分子結晶の候補を生成・緩和し、エネルギー順に並べます。``csp`` extra が必要です
（``pip install -e ".[csp]"``）。設定と出力は :doc:`csp_guide` を参照してください。

.. code-block:: bash

   molecular-modeling-csp [config]

**引数**:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - 引数
     - 説明
   * - ``config`` (任意)
     - CSP の設定ファイル（TOML）。省略時は ``input/csp.toml``

**出力**\ （``out_dir``、既定は ``csp_out/``）: ``candidates.csv``\ （順位表）、上位 50 件の
``NNN_<id>.cif``、``run.log``。``csp`` extra のパッケージが足りないときは、候補生成の前にその旨を出して止まります。

python -m molecular_modeling.csp.validate — 実験構造との比較
------------------------------------------------------------

CSP の実行結果の候補 CIF を実験 CIF と比べます（15 分子クラスターの RMSD）。
報告の読み方は :doc:`csp_guide` の「答え合わせ」を参照してください。

.. code-block:: bash

   python -m molecular_modeling.csp.validate <out_dir> <exp.cif> [<exp.cif> ...] [オプション]

**引数**:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - 引数
     - 説明
   * - ``out_dir``
     - ``molecular-modeling-csp`` の出力ディレクトリ（``candidates.csv`` と CIF）。
       存在しない、または候補がない場合はエラーで止まる
   * - ``exp.cif``
     - 実験 CIF ファイル（1 個以上。disorder のないもの。水素はなくてもよい）
   * - ``--n-top N``
     - 順位 ``N`` 以内の一致を成功とする（既定 20）
   * - ``--rmsd-cut X``
     - 一致とみなす RMSD\ :sub:`15` のしきい値（Å、既定 0.8）
   * - ``--flip A,B``
     - 別枠の flip 判定で同じ元素とみなす 2 元素（既定 ``O,S``。flip 判定は既定で有効）
   * - ``--no-flip``
     - flip 判定を省く（O や S を含まない分子で使う）

結果は画面と ``<out_dir>/validation.md`` に出ます。

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
