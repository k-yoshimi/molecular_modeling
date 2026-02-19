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

インストール手順
----------------

1. リポジトリのクローン
~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

   git clone <repository-url>
   cd molecular_modeling

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

4. 動作確認
~~~~~~~~~~~

テストを実行して正しくインストールされたことを確認します:

.. code-block:: bash

   pytest

全48テストがパスすれば、インストールは完了です。

Conda 環境の場合
~~~~~~~~~~~~~~~~

``environment.yml`` を使用して Conda 環境を構築できます:

.. code-block:: bash

   conda env create -f environment.yml
   conda activate molecular_modeling

詳細は ``CONDA_SETUP.md`` を参照してください。
