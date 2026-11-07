インストール
============

動作要件
--------

- Python 3.11 以上（3.13 で動作確認済み）
- pip または uv パッケージマネージャ

依存ライブラリ
~~~~~~~~~~~~~~

主要な依存ライブラリは以下の通りです（ ``pyproject.toml`` で管理）:

- **ASE** (3.25.0) — Atomic Simulation Environment
- **NumPy** (2.3.1) — 数値計算
- **SciPy** (1.16.0) — Voronoi 計算、最適化
- **matplotlib** (3.10.3) — 可視化
- **trimesh** (4.7.1) — メッシュ操作
- **manifold3d** (3.2.1) — マニフォールド再構成
- **Pillow** (11.3.0) — 画像処理

結晶構造予測（``csp`` extra。:doc:`csp_guide` を参照）で使う任意の依存ライブラリ:

- **PyXtal** (>= 1.1.5) — 対称性をもつ分子結晶のランダム生成
- **tblite** (>= 0.8.0) — 緩和に使う GFN2-xTB のエネルギーと力
- **pymatgen** (>= 2025.1.9) — CIF の読み書き。**spglib**\ （空間群の判定）も一緒に入る

機械学習ポテンシャル OrbMol-v2 で緩和する場合は、``csp-ml`` extra で次も入ります（Python 3.12 以上、約 1 GB）:

- **orb-models** (>= 0.7, < 0.8) — OrbMol-v2 のポテンシャル。**PyTorch** と **NVIDIA Warp** も一緒に入る

インストール手順
----------------

1. リポジトリのクローン
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone <repository-url>
   cd molecular_modeling   # リポジトリルート

2. 仮想環境の作成
~~~~~~~~~~~~~~~~~

.. code-block:: bash

   python -m venv .venv
   source .venv/bin/activate    # Linux / macOS
   .venv\Scripts\activate       # Windows

3. パッケージのインストール
~~~~~~~~~~~~~~~~~~~~~~~~~~~

pip を使用する場合:

.. code-block:: bash

   pip install -e .

uv を使用する場合:

.. code-block:: bash

   uv sync

結晶構造予測を使う場合は ``csp`` extra を加えます:

.. code-block:: bash

   pip install -e ".[csp]"     # uv なら: uv sync --extra csp
   pip install -e ".[csp-ml]"  # csp に OrbMol-v2 を加える（uv なら: uv sync --extra csp-ml）

4. 動作確認
~~~~~~~~~~~

テストを実行して正しくインストールされたことを確認します:

.. code-block:: bash

   pytest

全テストがパスすれば、インストールは完了です。

