クイックスタート
================

このガイドでは、最短手順でパイプラインを実行し、結果を確認する方法を説明します。

パイプラインの実行
------------------

仮想環境をアクティベートした状態で、以下を実行します:

.. code-block:: bash

   source .venv/bin/activate
   python run.py

実行が完了すると、以下のような出力が表示されます:

.. code-block:: text

   === Build system from TOML ===
   ...
   === Optimize Voronoi overlap ===
   ...
   === Write structure files ===
   Written: ./structure_out/qe.in
   Written: ./structure_out/POSCAR
   Written: ./structure_out/INCAR
   Written: ./structure_out/KPOINTS
   Written: ./structure_out/structure.cif

入力ファイル
------------

パイプラインは以下のファイルを読み込みます:

- ``input/input.toml`` — メインの設定ファイル（セルベクトル、分子リスト、繰り返し数など）
- ``data/<name>.toml`` — 分子定義ファイル（原子位置、格子ベクトルなど）

出力ファイル
------------

``structure_out/`` ディレクトリに以下のファイルが生成されます:

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - ファイル
     - 説明
   * - ``qe.in``
     - Quantum ESPRESSO 入力ファイル
   * - ``POSCAR``
     - VASP 構造ファイル
   * - ``INCAR``
     - VASP 計算設定ファイル
   * - ``KPOINTS``
     - VASP k点メッシュファイル
   * - ``structure.cif``
     - CIF 構造ファイル

``voronoi_pngs/`` ディレクトリに Voronoi セルの可視化画像が生成されます。

GUI の起動
----------

GUI を起動するには:

.. code-block:: bash

   python molecular_modeling_gui.py

既存の CIF ファイルを開く場合:

.. code-block:: bash

   python vis_cif.py structure.cif

GUI では以下の操作が可能です:

- 分子の追加・削除
- 元素の置換
- 構造ファイルのエクスポート
- 3D ビューの操作（回転・並進・ズーム）

詳細は :doc:`gui_guide` を参照してください。

次のステップ
------------

- :doc:`tutorial` — ベンゼン＋ナフタレンを使った詳細なチュートリアル
- :doc:`gui_guide` — GUI の全操作ガイド
- :doc:`cli_reference` — コマンドラインツールの使い方
- :doc:`toml_reference` — TOML 設定ファイルの全キー解説
- :doc:`architecture` — 内部アーキテクチャの説明
