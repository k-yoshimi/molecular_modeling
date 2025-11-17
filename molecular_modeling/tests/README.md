# ユニットテスト

このディレクトリには、`source`モジュールのユニットテストが含まれています。

## テストの実行

### すべてのテストを実行

```bash
pytest
```

### 特定のテストファイルを実行

```bash
pytest tests/test_build_molecule_groups.py
```

### 詳細な出力で実行

```bash
pytest -v
```

### カバレッジレポートを生成

```bash
pytest --cov=molecular_modeling.source --cov-report=html
```

## テストファイル一覧

- `test_build_molecule_groups.py`: 分子グループ構築のテスト
- `test_build_system_from_toml.py`: TOMLからのシステム構築のテスト
- `test_voronoi_geom.py`: Voronoi幾何計算のテスト
- `test_write_structure.py`: 構造ファイル出力のテスト
- `test_optimize_voronoi_overlap.py`: Voronoi最適化のテスト
- `test_render_voronoi_overlays.py`: Voronoi可視化のテスト
- `test_render_helpers.py`: レンダリングヘルパー関数のテスト

## フィクスチャ

`conftest.py`には以下の共通フィクスチャが定義されています：

- `simple_cell`: 単純な立方体セル
- `simple_atoms`: テスト用のAtomsオブジェクト
- `simple_polyhedron`: テスト用のPolyhedronオブジェクト

## テストの追加

新しいテストを追加する際は、以下の命名規則に従ってください：

- テストファイル: `test_<module_name>.py`
- テストクラス: `Test<ClassName>`
- テストメソッド: `test_<functionality_description>`

