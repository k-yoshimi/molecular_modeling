API リファレンス
================

このセクションでは、主要なモジュールと関数の API ドキュメントを自動生成しています。

モジュール構成
--------------

本ツールキットは以下のモジュール群で構成されています:

- **パイプライン**: ``main`` がエントリポイント。``build_system_from_toml`` → ``build_molecule_groups`` → ``optimize_voronoi_overlap`` → ``render_voronoi_overlays`` → ``write_structure`` の順で処理
- **幾何演算**: ``voronoi_geom`` が Voronoi セルの構築・体積計算・IoU 評価を担当。``Polyhedron`` クラスが中核
- **GUI**: ``gui.controller`` がビュー（ダイアログ）とモデル（ASE Atoms）を橋渡し。``gui.molecule_registry`` が ``data/*.toml`` の動的読み込みを管理
- **結晶構造予測**: ``csp.pipeline`` が ``csp.config`` の設定で生成・ふるい・緩和・順位付けを実行。``csp.validate`` が結果を実験 CIF と比較。``molecule_unwrap`` がセル境界で切れた分子をつなぎ直す（:doc:`../csp_guide` を参照）

プログラムから利用する場合の基本的な流れ:

.. code-block:: python

   from molecular_modeling.main import main

   # input/input.toml を読み込んでパイプライン全体を実行
   main()

メインパイプライン
------------------

main
~~~~

.. automodule:: molecular_modeling.main
   :members:
   :undoc-members:

システム構築
------------

build_system_from_toml
~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: molecular_modeling.build_system_from_toml
   :members:
   :undoc-members:

build_molecule_groups
~~~~~~~~~~~~~~~~~~~~~

.. automodule:: molecular_modeling.build_molecule_groups
   :members:
   :undoc-members:

Voronoi ジオメトリ
------------------

voronoi_geom
~~~~~~~~~~~~

.. automodule:: molecular_modeling.voronoi_geom
   :members:
   :undoc-members:
   :show-inheritance:

最適化
------

optimize_voronoi_overlap
~~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: molecular_modeling.optimize_voronoi_overlap
   :members:
   :undoc-members:

可視化
------

render_voronoi_overlays
~~~~~~~~~~~~~~~~~~~~~~~

.. automodule:: molecular_modeling.render_voronoi_overlays
   :members:
   :undoc-members:

出力
----

write_structure
~~~~~~~~~~~~~~~

.. automodule:: molecular_modeling.write_structure
   :members:
   :undoc-members:

分子と結晶構造予測
------------------

molecule_unwrap
~~~~~~~~~~~~~~~

.. automodule:: molecular_modeling.molecule_unwrap
   :members:
   :undoc-members:

csp.config
~~~~~~~~~~

.. automodule:: molecular_modeling.csp.config
   :members:
   :undoc-members:

csp.pipeline
~~~~~~~~~~~~

.. automodule:: molecular_modeling.csp.pipeline
   :members:
   :undoc-members:

csp.molecule
~~~~~~~~~~~~

.. automodule:: molecular_modeling.csp.molecule
   :members:
   :undoc-members:

csp.validate
~~~~~~~~~~~~

.. automodule:: molecular_modeling.csp.validate
   :members:
   :undoc-members:

GUI
---

controller
~~~~~~~~~~

.. automodule:: molecular_modeling.gui.controller
   :members:
   :undoc-members:

molecule_registry
~~~~~~~~~~~~~~~~~

.. automodule:: molecular_modeling.gui.molecule_registry
   :members:
   :undoc-members:
