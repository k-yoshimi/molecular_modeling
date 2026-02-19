Molecular Modeling Toolkit ドキュメント
==========================================

Molecular Modeling Toolkit は、ASE (Atomic Simulation Environment) をベースとした
分子モデリングツールキットです。TOML 形式の分子定義ファイルからスーパーセルを構築し、
Voronoi セル IoU 最大化による分子配向の最適化を行い、Quantum ESPRESSO / VASP / CIF
形式で出力します。

主な機能
--------

- **TOML ベースの分子定義**: 直感的な設定ファイルで分子構造を定義
- **スーパーセル構築**: 任意の繰り返し数でスーパーセルを生成
- **Voronoi 最適化**: Euler 角のグリッドサーチと Nelder-Mead 法による分子配向の最適化
- **4 種類の置換機能**: グローバル元素置換、分子置換、分子内原子置換、個別原子置換
- **マルチフォーマット出力**: QE, VASP, CIF 形式への一括エクスポート
- **GUI**: ASE-GUI を拡張したインタラクティブな操作環境

目次
----

.. toctree::
   :maxdepth: 2
   :caption: ユーザーガイド

   installation
   quickstart
   tutorial
   gui_guide

.. toctree::
   :maxdepth: 2
   :caption: リファレンス

   cli_reference
   toml_reference
   architecture
   api/index

他の言語
--------

- `English <../en/index.html>`_

索引
----

* :ref:`genindex`
* :ref:`modindex`
