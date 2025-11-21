# Conda環境セットアップ

このプロジェクト用のconda環境を作成・管理する方法です。

## 環境の作成

### 方法1: environment.ymlから作成（推奨）

```bash
cd molecular_modeling
conda env create -f environment.yml
```

### 方法2: 手動で作成

```bash
conda create -n molecular_modeling python=3.10
conda activate molecular_modeling
pip install -r requirements.txt
```

## 環境の有効化

```bash
conda activate molecular_modeling
```

## 環境の無効化

```bash
conda deactivate
```

## 環境の更新

依存関係が更新された場合：

```bash
conda env update -f environment.yml --prune
```

または、requirements.txtが更新された場合：

```bash
conda activate molecular_modeling
pip install -r requirements.txt --upgrade
```

## 環境の削除

```bash
conda env remove -n molecular_modeling
```

## 環境の一覧確認

```bash
conda env list
```

## 環境のエクスポート

現在の環境を`environment.yml`にエクスポートする場合：

```bash
conda env export > environment.yml
```

注意: エクスポートされたファイルには、システム固有のパスが含まれる場合があります。
共有用には`environment.yml`を手動で編集することを推奨します。

