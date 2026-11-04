やりたいこと別ガイド
====================

「こういうことがしたい」 → 「こうすればできる」の逆引きリファレンスです。
詳細は各リンク先を参照してください。

.. contents:: カテゴリ
   :local:
   :depth: 1

結晶構造を見たい・開きたい
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - やりたいこと
     - 操作
   * - CIF ファイルを開いて構造を見たい
     - GUI: **MoleculeTools → Load CIF...**
       またはコマンドラインで ``python vis_cif.py file.cif``
       （詳細は :doc:`cli_reference` を参照）
   * - QE / VASP / XYZ など他の形式を開きたい
     - GUI: **MoleculeTools → Load Structure...** (Ctrl+L)。
       ASE が対応する形式ならすべて読めます
   * - TOML プロジェクトを読み込みたい
     - GUI: **MoleculeTools → Load Project...** (Ctrl+T)。
       ``input.toml`` を選択すると分子グループ情報も復元されます
   * - 最近開いたファイルをすぐ開きたい
     - GUI: **Recent Files** メニューに最近の 10 件が表示されます

スーパーセルを作りたい
----------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - やりたいこと
     - 操作
   * - 単位胞を繰り返してスーパーセルにしたい
     - GUI: **MoleculeTools → Build Supercell...**
       で [a, b, c] 方向の繰り返し数を指定。
       詳細は :doc:`tutorial` のステップ 2 を参照
   * - セルパラメータ (a, b, c, α, β, γ) を確認・変更したい
     - GUI: **MoleculeTools → Cell Parameters...**
   * - TOML で一括設定したい
     - ``input.toml`` の ``supercell_repeat = [2, 2, 1]`` で指定。
       詳細は :doc:`toml_reference` を参照

分子を追加・配置したい
----------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - やりたいこと
     - 操作
   * - スーパーセルに分子を追加したい
     - GUI: **MoleculeTools → Add Molecule...** (Ctrl+J)。
       ``data/`` フォルダの TOML 分子定義から選択し、位置と配向を指定。
       詳細は :doc:`gui_guide` の Add Molecule ダイアログを参照
   * - 分子の位置を動かしたい
     - 分子中の原子を選択 → **Selection → Move Molecule/Atoms...**。
       X/Y/Z 方向の変位を ± ボタンで調整
   * - 分子の向きを回転したい
     - 分子中の原子を選択 → **Selection → Rotate Molecule/Atoms...**。
       トラックボールまたは Euler 角入力で回転
   * - 分子をマウスでドラッグ移動したい
     - **Selection → Toggle Drag Mode** で有効化。
       クリック＆ドラッグで直接移動
   * - 分子を削除したい
     - 分子中の原子を選択 → **Selection → Delete Molecule...**
   * - 選択した原子を新しい分子として保存したい
     - 原子を選択 → **MoleculeTools → Save Selection as Molecule...**。
       ``data/`` に TOML ファイルが作成されます

原子・元素を置換したい
----------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - やりたいこと
     - 操作
   * - 系全体で特定の元素を一括置換したい（例: 全 H → D）
     - GUI: **MoleculeTools → Global Element Substitution...**。
       または ``input.toml`` の ``[global_atom_substitutions]`` で設定
   * - 特定の分子だけ別の分子に丸ごと入れ替えたい
     - GUI: **MoleculeTools → Molecule Substitution...**。
       置換先の分子を選択すると、PCA で配向が維持されます
   * - 1 つの原子だけ別の元素に変えたい
     - GUI: **MoleculeTools → Individual Atom Substitution...**
   * - 複数の元素置換をまとめて実行したい（例: C→Si, H→F を同時に）
     - GUI: **MoleculeTools → Batch Element Substitution...**。
       複数のルールを一度に指定できます
   * - 特定の分子の全原子を新しい元素リストで再定義したい
     - GUI: **MoleculeTools → Molecule Atom Substitution...**

分子配向を最適化したい
----------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - やりたいこと
     - 操作
   * - Voronoi 最適化で分子配向を自動決定したい
     - GUI: **Pipeline → Run Optimize (Voronoi)** (Ctrl+P)。
       参照結晶との Voronoi セルの重なりを最大化します。
       アルゴリズムの詳細は :doc:`appendix_voronoi` を参照
   * - 最適化結果（IoU・回転角・スケール）を確認したい
     - GUI: **Pipeline → Show Optimization Results...** (Ctrl+H)。
       CSV エクスポートも可能です
   * - Voronoi セルの重なりを画像で確認したい
     - GUI: **Pipeline → Run Render (Voronoi PNGs)...**。
       青＝系、赤＝辞書、オレンジ＝交差部分
   * - 最適化から出力まで一括実行したい
     - GUI: **Pipeline → Run All (Optimize → Render → Export)...**
   * - コマンドラインで一括実行したい
     - ``python run.py``。
       詳細は :doc:`cli_reference` を参照

結晶構造を予測したい
--------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - やりたいこと
     - 操作
   * - 分子 1 個から結晶構造の候補を予測したい
     - ``molecular-modeling-csp input/csp.toml``\ （``pip install -e ".[csp]"`` が必要）。
       :doc:`csp_guide` を参照
   * - 既知の結晶から予測用の分子ファイルを作りたい
     - GUI: 分子 1 個を選び、**MoleculeTools → Save Selection as Molecule...**。
       :doc:`csp_guide` を参照
   * - 予測した候補を実験 CIF と比べたい
     - ``python -m molecular_modeling.csp.validate csp_out exp.cif``。
       :doc:`csp_guide` と :doc:`cli_reference` を参照

構造をファイルに出力したい
--------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - やりたいこと
     - 操作
   * - Quantum ESPRESSO の入力ファイルを作りたい
     - GUI: **Export → Quantum ESPRESSO...**
       または **Export → Export All (QE, VASP, CIF)...** (Ctrl+D)
   * - VASP (POSCAR/INCAR/KPOINTS) を作りたい
     - GUI: **Export → VASP...**
       または **Export → Export All...**
   * - CIF ファイルとして保存したい
     - GUI: **Export → CIF...**
       または **Export → Export All...**
   * - 3 形式すべてをまとめて出力したい
     - GUI: **Export → Export All (QE, VASP, CIF)...** (Ctrl+D)
   * - プロジェクト (TOML) として保存したい
     - GUI: **MoleculeTools → Save Project...** (Ctrl+K)。
       ``input.toml`` と ``data/*.toml`` が保存されます
   * - 3D ビューのスクリーンショットを保存したい
     - GUI: **View Settings → Save Screenshot...** (Ctrl+Shift+S)

表示を変えたい
--------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - やりたいこと
     - 操作
   * - 結合（ボンド）を表示したい
     - **View Settings → Show Bonds** (Ctrl+B) で切り替え。
       サイドパネルのチェックボックスでも操作可
   * - 水素原子を非表示にしたい
     - **View Settings → Toggle Hydrogen**
   * - 特定の元素だけ表示/非表示したい
     - **View Settings → Element Visibility...**
   * - 分子ごとに色分けしたい
     - **View Settings → Color by Molecule Group**
   * - 原子の色・サイズを変えたい
     - **View Settings → Atom Appearance...**
   * - 背景色を変えたい
     - ツールバーの **BG Color** ボタン（白/黒/ダーク/ライト/カスタム）
   * - ウィンドウサイズを変えたい
     - **View Settings → Small / Medium / Large / Wide / Custom Size...**
   * - サイドパネルを表示/非表示したい
     - **View Settings → Toggle Side Panel**

原子を選択したい
----------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - やりたいこと
     - 操作
   * - 全原子を選択したい
     - **Selection → Select All** (Ctrl+A)
   * - 選択を解除したい
     - **Selection → Deselect All** (Escape)
   * - 特定の元素だけ選択したい（例: 全 C）
     - **Selection → Select by Element...**
   * - 特定の分子グループだけ選択したい
     - **Selection → Select by Molecule...**
   * - ある原子の近傍を球状に選択したい
     - 原子を選択 → **Selection → Select by Sphere...**
   * - 選択した分子だけにフォーカスしたい
     - **Selection → Focus on Selection** (Ctrl+F)。
       元に戻すには Ctrl+Shift+F
   * - 原子間の距離を測りたい
     - 2 原子を選択するとステータスバーに距離が表示されます
   * - 結合角度を測りたい
     - 3 原子を選択するとステータスバーに角度が表示されます
   * - 二面角を測りたい
     - 4 原子を選択するとステータスバーに二面角が表示されます
   * - 測定履歴を見たい
     - **MoleculeTools → Measurement Log...**

操作を戻したい
--------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - やりたいこと
     - 操作
   * - 直前の操作を元に戻したい
     - **MoleculeTools → Undo** (Ctrl+Z)。最大 50 ステップ
   * - 戻した操作をやり直したい
     - **MoleculeTools → Redo** (Ctrl+Y)

TOML 設定ファイルの書き方を知りたい
------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 40 60

   * - やりたいこと
     - 操作
   * - input.toml の全パラメータを知りたい
     - :doc:`toml_reference` に全フィールドの説明があります
   * - 分子定義ファイル (data/\*.toml) の書き方を知りたい
     - :doc:`toml_reference` の「分子定義ファイル」セクションを参照
   * - 置換ルールの書き方を知りたい
     - :doc:`toml_reference` の「置換設定」セクションを参照
   * - 実際の使用例を見たい
     - :doc:`tutorial` にベンゼン＋ナフタレンの混合結晶の例があります
