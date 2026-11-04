# 元素置換分子の結晶構造予測 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 孤立分子1つと（任意で）置換元 CIF から、格子を含む結晶構造候補を生成し、GFN2-xTB のエネルギー順に並べて出力する CLI を作る。あわせて、既存 Voronoi パイプラインの2つのバグを直す。

**Architecture:** 新しいサブパッケージ `source/csp/` を作る。3つの生成器（PyXtal ランダム、Voronoi 向き、置換元由来）が `Candidate` を作り、ふるい分け（重なり、充填率、剛体 LJ）、tblite による GFN2-xTB 格子緩和、pymatgen による重複除去、出力の順に流す。既存の `voronoi_geom.py` は Voronoi セル計算として再利用する。既存の `run.py` パイプラインと GUI はそのまま残す。

**Tech Stack:** Python 3.11+（開発は 3.13）、ASE 3.29、PyXtal 1.1.5、tblite 0.8.0（GFN2-xTB。GFN-FF は無い）、pymatgen（StructureMatcher）、scipy、pytest。

**Spec:** `docs/superpowers/specs/2026-10-30-crystal-structure-prediction-design.md`

## Global Constraints

- ライセンスは LGPL-2.1-or-later。新しい依存は `[project.optional-dependencies] csp = ["pyxtal>=1.1.5", "tblite>=0.8.0", "pymatgen>=2025.1.9"]` に入れる。本体の `dependencies` には足さない。
- csp の依存を使うテストは、ファイルの先頭で `pytest.importorskip("pyxtal")` や `pytest.importorskip("tblite")` を呼ぶ。
- **原子座標を拡大・縮小しない。** 分子に許される操作は、剛体としての回転と平行移動だけ。例外は GFN2-xTB による緩和。
- 原子の並びは「分子ごとに連続」とし、ASE の `Atoms` に `mol_id`（int 配列）を持たせる。これは既存の規約と同じ。
- 分子の重心は、全体を通して**単純平均**（質量で重み付けしない）を使う。
- コードの規約は Ruff（E, F, W, I、行の長さ 120）、mypy（`source/csp/` は対象、`source/gui/` は対象外）、numpy 形式の docstring。これは既存コードに合わせる。
- 時間の目安は既定で `time_budget_min = 60`。
- 既定の空間群は `[14, 2, 19, 15, 4, 61]`（P2₁/c、P-1、P2₁2₁2₁、C2/c、P2₁、Pbca）。
- しきい値の既定値:
  - 重なり判定の係数 `overlap_factor = 0.7`（ファンデルワールス半径の和に掛ける）
  - 充填率 `packing_range = [0.6, 0.8]`
  - 緩和後の分子 RMSD 上限 `max_mol_rmsd = 0.5` Å（水素を除く）
  - 答え合わせ: 上位 `n_top = 20` 件、`rmsd_cut = 0.8` Å、クラスターは 15 分子
- 作業ブランチ: `main` から `feature/csp` を切る。今の `pr-2` にある未コミットの変更（CLAUDE.md、docs/*/gui_guide.rst、uv.lock など）は持ち込まない。

## Review Focus

1. **GFN2-xTB の計算が収束しない、または例外を投げる候補。** 1件の失敗でパイプライン全体が止まってはいけない。その候補を `status="failed"` として記録し、残りの処理を続ける。→ Task 8 にテスト `test_relax_many_survives_failure` を追加する。
2. **分子 TOML が周期境界で切れたまま保存されている（加藤さんの既存ファイル）。** `load_molecule` は読み込みの時点でつながりを確かめ、切れていれば分かりやすいエラーを出す。黙って断片のまま使ってはいけない。→ Task 3 にテスト `test_load_molecule_rejects_fragmented` を追加する。
3. **置換元 CIF に部分占有（disorder）が残っている。** 対応しないので、どのサイトかを示すエラーを出す。→ Task 5 にテスト `test_seed_rejects_partial_occupancy` を追加する。
4. **ある空間群で PyXtal の生成が失敗し続ける**（分子が大きく、ワイコフ位置が合わないなど）。その空間群は数を減らしてログに残し、ほかの空間群は続ける。→ Task 4 にテスト `test_generate_random_skips_impossible_group` を追加する。
5. **すべての候補がふるい分けで落ちて、0 件になる。** 空の CSV と、「各段で何件残ったか」のログを出して、正常に終了する。→ Task 9 にテスト `test_run_csp_handles_zero_survivors` を追加する。

---

## ファイル構成

| ファイル | 役割 |
|---|---|
| `source/molecule_unwrap.py`（新規） | 周期境界で切れた分子をつなぎ直す。GUI の保存処理と置換元の読み込みで共用する |
| `source/optimize_voronoi_overlap.py`（修正） | 原子座標への倍率 `s` を外す |
| `source/gui/save_molecule_dialog.py`（修正） | `extract_molecule` を使って保存する |
| `source/csp/__init__.py` | 空のファイル |
| `source/csp/config.py` | `CSPConfig` と `load_csp_config` |
| `source/csp/molecule.py` | 分子の読み込み、原子の対応付けと重ね合わせ（Kabsch）、RMSD |
| `source/csp/candidate.py` | `Candidate` データクラス |
| `source/csp/gen_random.py` | PyXtal によるランダム生成 |
| `source/csp/gen_seed.py` | 置換元由来の生成 |
| `source/csp/gen_voronoi.py` | Voronoi で向きを決める生成 |
| `source/csp/screen.py` | 重なり、充填率、剛体 LJ によるふるい分け |
| `source/csp/relax.py` | 計算手段の作成（差し替え可能）、GFN2-xTB 緩和、並列実行 |
| `source/csp/dedupe.py` | 重複除去 |
| `source/csp/pipeline.py` | `run_csp`（全体の流れと出力） |
| `source/csp/validate.py` | 答え合わせ（RMSD₁₅）と、その CLI |
| `source/csp/cli.py` | `molecular-modeling-csp` の入口 |
| `input/csp.toml` | 入力例（ベンゼン） |
| `tests/csp/…` | 各モジュールのテスト |
| `docs/ja/csp_guide.rst`、`docs/en/csp_guide.rst` | 使い方 |

---

### Task 1: Voronoi 最適化で原子を拡大・縮小しないようにする

**Files:**
- Modify: `source/optimize_voronoi_overlap.py`（docstring の Notes の 3〜5 と、183 行目付近の `pos[idxs] = p_ref + s * (...)`）
- Test: `tests/test_optimize_voronoi_overlap.py`

**Interfaces:**
- Produces: `optimize_voronoi_overlap` のシグネチャは変えない。`system.info["optimized_scales"]` は Voronoi セルの比較に使った倍率として残す（ログの互換性のため）。

- [ ] **Step 0: ブランチを作る**

```bash
git switch main && git switch -c feature/csp
git add docs/superpowers/specs/2026-10-30-crystal-structure-prediction-design.md docs/superpowers/plans/2026-10-30-crystal-structure-prediction.md
git commit -m "docs: add design and plan for crystal structure prediction"
```

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_optimize_preserves_intramolecular_distances(self, simple_atoms):
    # Build a system whose cell volume differs from the dictionary lattice so that s != 1.
    # (Reuse the setup from test_optimize_voronoi_overlap_basic but give dict_meta a lattice 1.3x larger.)
    before = [simple_atoms.get_distances(g[0], g[1:]) for g in groups if len(g) > 1]
    result = optimize_voronoi_overlap(simple_atoms, groups=groups, grid_num=(3, 3, 3), rng=1)
    after = [result.get_distances(g[0], g[1:]) for g in groups if len(g) > 1]
    assert not np.allclose(result.info["optimized_scales"], 1.0)
    for b, a in zip(before, after):
        np.testing.assert_allclose(a, b, atol=1e-8)
```

テストでは、各サイトが2原子以上になるようにする。conftest の `simple_atoms` が1サイト1原子なら、テストの中で2原子の分子を組み立てる。

- [ ] **Step 2: 失敗を確認する**

Run: `pytest tests/test_optimize_voronoi_overlap.py::TestOptimizeVoronoiOverlap::test_optimize_preserves_intramolecular_distances -v`（クラス名は実際のファイルに合わせる）
Expected: FAIL（距離が `s` 倍になっている）

- [ ] **Step 3: 実装する**

`pos[idxs] = p_ref + rot.apply(sub - p_ref)` に変える。docstring の Notes の手順 5 を「剛体として回転だけを適用する。倍率 s は Voronoi セルの比較にだけ使う」に直す。

- [ ] **Step 4: テストを通す**

Run: `pytest tests/test_optimize_voronoi_overlap.py -v`
Expected: 全件 PASS

- [ ] **Step 5: コミットする**

```bash
git add source/optimize_voronoi_overlap.py tests/test_optimize_voronoi_overlap.py
git commit -m "fix: keep molecules rigid in Voronoi orientation optimization"
```

---

### Task 2: 周期境界で切れた分子をつなぎ直して保存する

**Files:**
- Create: `source/molecule_unwrap.py`
- Modify: `source/gui/save_molecule_dialog.py`（`apply` の 336〜352 行目）
- Test: `tests/test_molecule_unwrap.py`

**Interfaces:**
- Produces:
  - `unwrap_molecule(positions: np.ndarray, symbols: list[str], cell: np.ndarray, pbc: bool | Sequence[bool] = True) -> np.ndarray`。戻り値は (N,3) の直交座標。最初の原子を固定し、結合グラフを幅優先でたどって、各原子を最小像規約で隣の原子の近くに置き直す。
  - `extract_molecule(atoms: ase.Atoms, indices: Sequence[int]) -> tuple[list[str], np.ndarray, np.ndarray]`。戻り値は（記号、重心を原点にした座標、重心）。重心は単純平均。
  - `is_connected(positions: np.ndarray, symbols: list[str]) -> bool`。周期境界なしで、結合グラフが1つにつながっているかを返す。

結合判定には `ase.neighborlist.natural_cutoffs(..., mult=1.2)` を使う。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_unwrap_rejoins_molecule_split_across_boundary():
    cell = np.diag([10.0, 10.0, 10.0])
    # C-C 1.5 Å bond straddling x=0: one atom at x=9.4, the other at x=0.9 (wrapped)
    pos = np.array([[9.4, 5.0, 5.0], [0.9, 5.0, 5.0]])
    out = unwrap_molecule(pos, ["C", "C"], cell)
    assert np.isclose(np.linalg.norm(out[1] - out[0]), 1.5)
    assert is_connected(out, ["C", "C"])
    assert not is_connected(pos, ["C", "C"])

def test_extract_molecule_uses_plain_mean_and_unwraps():
    atoms = Atoms("CO", positions=[[9.4, 5, 5], [0.7, 5, 5]], cell=np.diag([10.0] * 3), pbc=True)
    syms, centered, com = extract_molecule(atoms, [0, 1])
    assert syms == ["C", "O"]
    np.testing.assert_allclose(centered.mean(axis=0), 0.0, atol=1e-12)
    np.testing.assert_allclose(com, [(9.4 + 10.7) / 2, 5, 5])
```

- [ ] **Step 2: 失敗を確認する**

Run: `pytest tests/test_molecule_unwrap.py -v`
Expected: FAIL（`ModuleNotFoundError: molecular_modeling.molecule_unwrap`）

- [ ] **Step 3: `source/molecule_unwrap.py` を実装し、GUI の `apply` を `extract_molecule` を使う形に置き換える**

`apply` は `extract_molecule` を呼び、返ってきた値を `_format_toml(positions=centered, center_of_mass=com[None, :])` に渡す。`sub_atoms.get_center_of_mass()` の呼び出しは消す。

- [ ] **Step 4: テストを通す**

Run: `pytest tests/test_molecule_unwrap.py tests/test_gui_features.py -v`
Expected: PASS

- [ ] **Step 5: コミットする**

```bash
git add source/molecule_unwrap.py source/gui/save_molecule_dialog.py tests/test_molecule_unwrap.py
git commit -m "fix: rejoin molecules split by periodic boundaries when saving"
```

---

### Task 3: csp の土台（依存、設定、分子、Candidate、重ね合わせ）

**Files:**
- Modify: `pyproject.toml`（`csp` extra、`[project.scripts] molecular-modeling-csp = "molecular_modeling.csp.cli:main"`）
- Create: `source/csp/__init__.py`、`source/csp/config.py`、`source/csp/molecule.py`、`source/csp/candidate.py`、`input/csp.toml`
- Test: `tests/csp/__init__.py`、`tests/csp/conftest.py`（`default_cfg() -> CSPConfig`。ベンゼン、`space_groups=[2]`、`n_candidates_per_sg=10`、`n_procs=1`、`relax_steps=5`、`out_dir=Path("csp_out_test")` で `CSPConfig` を直接組み立てて返す普通の関数。後の Task のテストが import する）、`tests/csp/test_config.py`、`tests/csp/test_molecule.py`

**Interfaces:**
- Consumes: `is_connected`（Task 2）
- Produces:
  - `@dataclass CSPConfig`。フィールド: `molecule: Path`、`mode: Literal["blind", "informed"]`、`space_groups: list[int]`、`z: dict[int, int] | None`（空間群 → Z。`None` のときは一般位置の多重度 × `z_prime`）、`z_prime: int = 1`、`n_candidates_per_sg: int = 500`、`generators: dict[str, float]`（既定 `{"random": 0.6, "voronoi": 0.4}`。`seed_structures` があると `"seed"` を足す）、`seed_structures: list[SeedSpec]`、`packing_range: tuple[float, float] = (0.6, 0.8)`、`overlap_factor: float = 0.7`、`n_prescreen_keep: int = 150`、`relax_steps: int = 150`、`relax_fmax: float = 0.05`、`n_procs: int = os.cpu_count()`、`max_mol_rmsd: float = 0.5`、`time_budget_min: float = 60`、`rng: int = 1`、`out_dir: Path = Path("csp_out")`。
  - `@dataclass SeedSpec(cif: Path, substitutions: list[tuple[str, str, int]], volume_scale: float | Literal["auto"] = "auto")`。`substitutions` の各要素は（元の元素、新しい元素、1分子あたりの置換数）。
  - `load_csp_config(path: str | Path) -> CSPConfig`。`mode="blind"` のときは `space_groups` を既定リストで上書きし、`z` は `None` にする。`informed` のときに `space_groups` が無ければ `ValueError`。
  - `load_molecule(path: str | Path) -> ase.Atoms`。`.toml`（既存の `data/*.toml` 形式。`use_fractional` があれば `lattice_vectors` で変換）と `.xyz` を読む。重心を原点に移す。つながっていなければ `ValueError("molecule ... is fragmented ...")`。
  - `match_and_align(mobile: ase.Atoms, target: ase.Atoms, heavy_only: bool = True) -> AlignResult`。`AlignResult(rotation: np.ndarray (3,3), translation: np.ndarray (3,), perm: np.ndarray, rmsd: float)`。`perm[i]` は target の i 番目の原子に対応する mobile の原子の番号。処理は次の通り。(1) 主軸をそろえる初期値を4通り（符号の組み合わせのうち det=+1 のもの）作る。(2) 元素ごとに `scipy.optimize.linear_sum_assignment` で対応付ける。(3) Kabsch で重ねる。(2) と (3) を収束まで（最大 20 回）くり返し、RMSD が最小のものを返す。元素の内訳が違う場合（O→S 置換の前後）は、`element_map: dict[str, str] | None` 引数で同じ元素として扱う組を与える。
  - `molecule_rmsd(ref: ase.Atoms, mol: ase.Atoms) -> float`。`match_and_align(...).rmsd` の別名。
  - `@dataclass Candidate(atoms: ase.Atoms, origin: str, space_group: int, meta: dict = field(default_factory=dict), energy_per_mol: float | None = None, density: float | None = None, status: str = "generated")`。`origin` は `"random"`、`"voronoi"`、`"seed:<cifの名前>"` のいずれか。`atoms.arrays["mol_id"]` が必須。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_blind_mode_uses_default_space_groups(tmp_path):
    p = tmp_path / "csp.toml"
    p.write_text('molecule = "data/benzene.toml"\nmode = "blind"\nspace_groups = [1]\n')
    cfg = load_csp_config(p)
    assert cfg.space_groups == [14, 2, 19, 15, 4, 61]
    assert cfg.z is None

def test_informed_mode_requires_space_groups(tmp_path):
    p = tmp_path / "csp.toml"
    p.write_text('molecule = "data/benzene.toml"\nmode = "informed"\n')
    with pytest.raises(ValueError):
        load_csp_config(p)

def test_load_molecule_centers_benzene():
    mol = load_molecule("data/benzene.toml")
    assert len(mol) == 12
    np.testing.assert_allclose(mol.get_positions().mean(axis=0), 0.0, atol=1e-10)

def test_load_molecule_rejects_fragmented(tmp_path):  # Review Focus 2
    # two C atoms 6 Å apart -> not bonded
    (tmp_path / "frag.xyz").write_text("2\n\nC 0 0 0\nC 6 0 0\n")
    with pytest.raises(ValueError, match="fragmented"):
        load_molecule(tmp_path / "frag.xyz")

def test_match_and_align_recovers_rotation_and_permutation():
    mol = load_molecule("data/benzene.toml")
    rot = Rotation.from_euler("xyz", [30, 60, 10], degrees=True)
    moved = mol.copy(); moved.positions = rot.apply(mol.positions) + [1, 2, 3]
    moved = moved[np.random.default_rng(0).permutation(len(moved))]
    res = match_and_align(moved, mol)
    assert res.rmsd < 1e-6
```

- [ ] **Step 2: 失敗を確認する**

Run: `pytest tests/csp -v`
Expected: FAIL（`ModuleNotFoundError`）

- [ ] **Step 3: 実装し、`pip install -e ".[dev,csp]"` を実行する**

`input/csp.toml` には、ベンゼン、blind モード、`n_candidates_per_sg = 50` の例を書き、各キーに1行ずつコメントを付ける。

- [ ] **Step 4: テストを通し、パッケージとして読み込めることを確認する**

Run: `pytest tests/csp -v && python -c "import molecular_modeling.csp.config"`
Expected: PASS、エラーなし。`csp` サブパッケージが読み込めない場合は、`pyproject.toml` の `[tool.setuptools]` に `packages = ["molecular_modeling", "molecular_modeling.gui", "molecular_modeling.csp"]` を足す。

- [ ] **Step 5: コミットする**

```bash
git add pyproject.toml source/csp input/csp.toml tests/csp
git commit -m "feat: add CSP config, molecule loading and alignment"
```

---

### Task 4: PyXtal によるランダム生成

**Files:**
- Create: `source/csp/gen_random.py`
- Test: `tests/csp/test_gen_random.py`

**Interfaces:**
- Consumes: `Candidate`、`CSPConfig`、`load_molecule`
- Produces:
  - `to_pyxtal_molecule(mol: ase.Atoms) -> pyxtal.molecule.pyxtal_molecule`（`pymatgen.core.Molecule` 経由）
  - `generate_random(mol: ase.Atoms, space_group: int, z: int, n: int, rng: np.random.Generator, max_attempts_factor: int = 5) -> list[Candidate]`。返す候補には `meta["pyxtal"]`（`pyxtal` オブジェクト。Task 6 で向きを書き換えるのに使う）を入れる。`atoms` は `crystal.to_ase(resort=False)` で作る。`mol_id` は、原子番号を分子の原子数で割った商で付ける。`n * max_attempts_factor` 回試しても `n` 件に届かなければ、作れた分だけ返し、`logging.warning` で空間群と件数を出す。
  - `default_z(space_group: int, z_prime: int) -> int`（一般位置の多重度 × `z_prime`。`pyxtal.symmetry.Group(sg)[0].multiplicity` を使う）

PyXtal の呼び出し方は `c = pyxtal(molecular=True); c.from_random(3, sg, [pmol], [z], seed=int)`。生成のあと `c.valid` を確かめる（2026-10-30 に動作を確認済み）。

- [ ] **Step 1: 失敗するテストを書く**

```python
pytest.importorskip("pyxtal")

def test_generate_random_benzene_p21c():
    mol = load_molecule("data/benzene.toml")
    cands = generate_random(mol, 14, 4, n=3, rng=np.random.default_rng(1))
    assert len(cands) == 3
    for c in cands:
        assert c.origin == "random" and c.space_group == 14
        assert len(c.atoms) == 48
        assert sorted(set(c.atoms.arrays["mol_id"])) == [0, 1, 2, 3]
        assert molecule_rmsd(mol, c.atoms[c.atoms.arrays["mol_id"] == 0]) < 1e-3  # rigid

def test_generate_random_skips_impossible_group(caplog):  # Review Focus 4
    mol = load_molecule("data/benzene.toml")
    cands = generate_random(mol, 230, 1, n=2, rng=np.random.default_rng(1), max_attempts_factor=1)
    assert cands == []
    assert "230" in caplog.text

def test_default_z():
    assert default_z(14, 1) == 4 and default_z(2, 2) == 4
```

- [ ] **Step 2: 失敗を確認する**

Run: `pytest tests/csp/test_gen_random.py -v`
Expected: FAIL

- [ ] **Step 3: 実装する。** PyXtal が例外を投げた場合も「失敗した1回」として数える。

- [ ] **Step 4: テストを通す**

Run: `pytest tests/csp/test_gen_random.py -v`
Expected: PASS

- [ ] **Step 5: コミットする**

```bash
git add source/csp/gen_random.py tests/csp/test_gen_random.py
git commit -m "feat: generate random molecular crystal candidates with PyXtal"
```

---

### Task 5: 置換元由来の生成

**Files:**
- Create: `source/csp/gen_seed.py`
- Test: `tests/csp/test_gen_seed.py`、`tests/csp/data/seed_benzene_p21c.cif`（Task 4 の生成器で作ったベンゼンの P2₁/c 結晶を、原子をセル内に折り返して CIF に書いたもの）

**Interfaces:**
- Consumes: `unwrap_molecule`、`match_and_align`、`Candidate`、`SeedSpec`
- Produces:
  - `read_seed_molecules(cif: Path) -> tuple[np.ndarray, list[ase.Atoms]]`。戻り値は（セル、つなぎ直した分子のリスト）。ASE で読み、占有率が 1 未満のサイトがあれば `ValueError` を出し、サイトのラベルを含める。周期境界ありの結合グラフで分子に分け、それぞれ `unwrap_molecule` でつなぎ直す。
  - `generate_from_seed(mol: ase.Atoms, seed: SeedSpec) -> list[Candidate]`。手順は次の通り。(1) 置換元の分子0について、`substitutions` を満たす原子の選び方を列挙する（たとえば O が2つで1つを S にするなら2通り）。(2) 各選び方について、分子0の置換後の構造に入力分子を `match_and_align`（`element_map` なし）で重ねる。(3) 分子 k ≥ 1 には、分子0から分子 k への剛体変換（置換元の分子同士を `match_and_align` して得る）を、重ねた入力分子に当てはめる。これで全分子の向きがそろう。(4) 格子は置換元のまま。`volume_scale="auto"` のときは、(入力分子の体積 / 置換元の分子の体積) の比で格子を等方的にスケールし、重心だけを移動する（分子は剛体のまま）。返すのは、置換の選び方の数 × 1 件。`origin = f"seed:{cif.stem}"`、`space_group = 1`、`meta["substitution_choice"]` に選んだ原子の番号を入れる。
  - 分子の体積は、ファンデルワールス半径の球の和を、格子点上のモンテカルロ（1 Å³ に 50 点程度）で見積もる。`screen.py` でも使うので、`molecule_volume(mol: ase.Atoms) -> float` として `source/csp/molecule.py` に置く。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_seed_rejoins_and_keeps_rigid_molecules():
    mol = load_molecule("data/benzene.toml")
    cands = generate_from_seed(mol, SeedSpec(cif=DATA / "seed_benzene_p21c.cif", substitutions=[], volume_scale=1.0))
    assert len(cands) == 1
    atoms = cands[0].atoms
    for k in set(atoms.arrays["mol_id"]):
        assert molecule_rmsd(mol, atoms[atoms.arrays["mol_id"] == k]) < 1e-3

def test_seed_substitution_enumerates_uniform_choices():
    # substitute 1 of benzene's 6 H by F: 6 choices, symmetric copies uniform across molecules
    mol = load_molecule("data/benzene.toml"); mol_f = mol.copy(); mol_f.symbols[6] = "F"
    cands = generate_from_seed(mol_f, SeedSpec(cif=DATA / "seed_benzene_p21c.cif", substitutions=[("H", "F", 1)], volume_scale=1.0))
    assert len(cands) == 6
    assert all(c.atoms.get_chemical_symbols().count("F") == 4 for c in cands)

def test_seed_rejects_partial_occupancy(tmp_path):  # Review Focus 3
    cif = write_seed_cif(tmp_path / "dis.cif", occupancy_of_first_site=0.7)
    with pytest.raises(ValueError, match="occupancy"):
        read_seed_molecules(cif)
```

`write_seed_cif(path, occupancy_of_first_site=1.0) -> Path` はテストファイルの中のヘルパー。Task 4 の `generate_random(benzene, 14, 4, n=1, rng=default_rng(0))` の結晶を `wrap()` してから、`_atom_site_occupancy` 列付きの CIF を手で書き出す。最初のサイトの占有率だけを引数の値にする。`tests/csp/data/seed_benzene_p21c.cif` も、このヘルパーで一度だけ作ってコミットする。

- [ ] **Step 2: 失敗を確認する**

Run: `pytest tests/csp/test_gen_seed.py -v`
Expected: FAIL

- [ ] **Step 3: 実装する**

- [ ] **Step 4: テストを通す**

Run: `pytest tests/csp/test_gen_seed.py -v`
Expected: PASS

- [ ] **Step 5: コミットする**

```bash
git add source/csp/gen_seed.py source/csp/molecule.py tests/csp/test_gen_seed.py tests/csp/data
git commit -m "feat: derive candidates from substituted seed crystals"
```

---

### Task 6: Voronoi で向きを決める生成

**Files:**
- Create: `source/csp/gen_voronoi.py`
- Test: `tests/csp/test_gen_voronoi.py`

**Interfaces:**
- Consumes: `voronoi_geom._build_system_seeds_cart`、`voronoi_geom._voronoi_cell_vertices_halfspaces_cart`、Task 4 の候補の `meta["pyxtal"]`、Task 5 の候補
- Produces:
  - `molecule_shape_points(mol: ase.Atoms, n_sphere: int = 26) -> np.ndarray`。各原子のまわりに、ファンデルワールス半径の球面上の点を `n_sphere` 個置いた点群を返す（重心が原点）。分子の形はこの点群の凸包で表す。
  - `site_voronoi_vertices(cell: np.ndarray, coms: np.ndarray, index: int, rng: int = 1) -> np.ndarray`。単位胞で切り取らない Voronoi セルの頂点を返す（重心が原点）。
  - `convex_iou(a_pts: np.ndarray, b_pts: np.ndarray) -> float`。2つの点群の凸包の IoU。共通部分は、両方の `ConvexHull.equations` をまとめて `HalfspaceIntersection` で求める（内点は原点。両方とも原点を含む）。体積は `ConvexHull(...).volume`。
  - `best_orientation(shape_pts: np.ndarray, cell_vertices: np.ndarray, grid_num: tuple[int, int, int] = (9, 9, 9)) -> tuple[Rotation, float]`。形の点群を `(V_cell / V_shape)^(1/3)` 倍して（比較のためだけ）、Euler 角のグリッド探索のあと Nelder–Mead で IoU を最大にする。既存の `optimize_voronoi_overlap` と同じ流れ。
  - `orient_candidate(cand: Candidate) -> list[Candidate]`。
    - `meta["pyxtal"]` がある場合：非対称単位の分子（`mol_id == 0`）について `best_orientation` で Δ を求める。`site.orientation.matrix = Δ.as_matrix() @ site.orientation.matrix` として、`to_ase(resort=False)` で作り直す。
    - ない場合（置換元由来）：分子ごとに Δ を求めて回転させる。
    - 返すのは2件。1件目は Δ。2件目は「Δ のあと、分子の最大慣性主軸（最も細長い方向）のまわりに 180° 回したもの」。`origin = "voronoi"`、`meta["voronoi_iou"]` に IoU を入れる。
  - `generate_voronoi(skeletons: list[Candidate]) -> list[Candidate]`。`orient_candidate` を各骨格に当てはめて、結果をつなげる。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_convex_iou_identity_and_disjoint():
    cube = np.array(list(itertools.product([-1, 1], repeat=3)), float)
    assert convex_iou(cube, cube) == pytest.approx(1.0)
    assert convex_iou(cube, cube * 0.5) == pytest.approx(1 / 8)

def test_best_orientation_aligns_elongated_shape_with_elongated_cell():
    # a rod along x in an orthorhombic cell elongated along z -> rod should end up along z
    rod = np.array([[x, y, z] for x in (-3, 3) for y in (-0.5, 0.5) for z in (-0.5, 0.5)], float)
    cell = np.diag([4.0, 4.0, 12.0]); coms = np.zeros((1, 3))
    rot, iou = best_orientation(rod, site_voronoi_vertices(cell, coms, 0))
    assert abs(rot.apply([1, 0, 0])[2]) > 0.95

def test_orient_candidate_keeps_symmetry_and_rigidity():
    mol = load_molecule("data/benzene.toml")
    skel = generate_random(mol, 14, 4, n=1, rng=np.random.default_rng(2))[0]
    out = orient_candidate(skel)
    assert len(out) == 2 and all(c.origin == "voronoi" for c in out)
    for c in out:
        for k in range(4):
            assert molecule_rmsd(mol, c.atoms[c.atoms.arrays["mol_id"] == k]) < 1e-3
        assert spglib_number(c.atoms) == 14  # helper in test using spglib.get_symmetry_dataset(symprec=1e-3)
```

- [ ] **Step 2: 失敗を確認する**

Run: `pytest tests/csp/test_gen_voronoi.py -v`
Expected: FAIL

- [ ] **Step 3: 実装する**

- [ ] **Step 4: テストを通す**

Run: `pytest tests/csp/test_gen_voronoi.py -v`
Expected: PASS

- [ ] **Step 5: コミットする**

```bash
git add source/csp/gen_voronoi.py tests/csp/test_gen_voronoi.py
git commit -m "feat: orient molecules by matching molecular shape to Voronoi cells"
```

---

### Task 7: ふるい分け（重なり、充填率、剛体 LJ）

**Files:**
- Create: `source/csp/screen.py`
- Test: `tests/csp/test_screen.py`

**Interfaces:**
- Consumes: `Candidate`、`molecule_volume`
- Produces:
  - `has_intermolecular_overlap(atoms: ase.Atoms, factor: float = 0.7) -> bool`。分子が違う原子の組で、距離がファンデルワールス半径の和 × `factor` 未満のものがあれば True。周期境界を考慮する（`ase.neighborlist.neighbor_list`）。半径は `ase.data.vdw_radii`。
  - `packing_fraction(atoms: ase.Atoms, mol_volume: float) -> float`。`n_mol * mol_volume / cell_volume`。
  - `rigid_lj_energy(atoms: ase.Atoms, cutoff: float = 10.0) -> float`。分子間の原子の組だけを足した、1分子あたりの LJ エネルギー（eV）。パラメータは UFF の (x_i, D_i)。H、C、N、O、F、S、Cl、Se、Br、I の値を表として持ち、組み合わせ則は幾何平均。表にない元素は `KeyError` を出し、元素名を含める。
  - `prescreen(cands: list[Candidate], mol_volume: float, cfg: CSPConfig) -> tuple[list[Candidate], dict[str, int]]`。重なり → 充填率 → LJ の順でふるい分け、LJ の低い順に `cfg.n_prescreen_keep` 件を返す。2つ目の戻り値は、各段で残った件数。各候補の `meta["lj"]` と `density` を埋める。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_overlap_detected_for_squeezed_cell():
    c = generate_random(load_molecule("data/benzene.toml"), 2, 2, n=1, rng=np.random.default_rng(0))[0]
    squeezed = c.atoms.copy(); squeezed.set_cell(squeezed.cell * 0.7, scale_atoms=False)
    # scale only COMs is not needed here: shrinking the cell without moving atoms brings periodic images closer
    assert has_intermolecular_overlap(squeezed)

def test_lj_prefers_reasonable_density():
    mol = load_molecule("data/benzene.toml")
    c = generate_random(mol, 2, 2, n=1, rng=np.random.default_rng(0))[0]
    e0 = rigid_lj_energy(c.atoms)
    expanded = expand_coms(c.atoms, 1.5)  # test helper: scale cell by 1.5 and move each molecule rigidly by its COM
    assert rigid_lj_energy(expanded) > e0  # a sparse crystal binds less

def test_prescreen_counts_and_keep():
    mol = load_molecule("data/benzene.toml")
    cands = generate_random(mol, 2, 2, n=10, rng=np.random.default_rng(0))
    kept, counts = prescreen(cands, molecule_volume(mol), replace(default_cfg(), n_prescreen_keep=3))
    assert len(kept) <= 3
    assert set(counts) == {"generated", "no_overlap", "packing", "kept"} and counts["generated"] == 10
    assert all(k.meta["lj"] <= c.meta.get("lj", np.inf) for k in kept for c in kept[kept.index(k):])  # sorted ascending
```

`expand_coms(atoms, factor) -> Atoms` はテストファイルの中のヘルパー。格子を `factor` 倍し、各分子を重心ごと（分率座標を保って）剛体として動かす。

- [ ] **Step 2: 失敗を確認する**

Run: `pytest tests/csp/test_screen.py -v`
Expected: FAIL

- [ ] **Step 3: 実装する**

- [ ] **Step 4: テストを通す**

Run: `pytest tests/csp/test_screen.py -v`
Expected: PASS

- [ ] **Step 5: コミットする**

```bash
git add source/csp/screen.py tests/csp/test_screen.py
git commit -m "feat: add overlap, packing and rigid LJ prescreening"
```

---

### Task 8: GFN2-xTB による緩和（並列、差し替え可能）

**Files:**
- Create: `source/csp/relax.py`
- Test: `tests/csp/test_relax.py`

**Interfaces:**
- Consumes: `Candidate`、`molecule_rmsd`、`is_connected`
- Produces:
  - `make_calculator(method: str = "GFN2-xTB") -> ase.calculators.calculator.Calculator`。今の対応は `"GFN2-xTB"` と `"GFN1-xTB"`（どちらも `tblite.ase.TBLite(method=..., verbosity=0)`）。テスト用に `"EMT"` も受け付ける。それ以外は `ValueError`。
  - `relax_candidate(cand: Candidate, ref_mol: ase.Atoms, method: str, fmax: float, steps: int, max_mol_rmsd: float) -> Candidate`。`FrechetCellFilter` と `FIRE` で緩和する。終わったら `energy_per_mol`、`density`（g/cm³）、`meta["max_mol_rmsd"]`、`meta["converged"]` を埋める。緩和後の分子のどれかが `max_mol_rmsd` を超えるか、つながりが切れていたら `status="rejected_deformed"`。例外が出たら `status="failed"` にして `meta["error"]` に例外の文字列を入れる。それ以外は `status="relaxed"`。例外は外に投げない。
  - `relax_many(cands: list[Candidate], ref_mol: ase.Atoms, cfg: CSPConfig, method: str = "GFN2-xTB") -> list[Candidate]`。`concurrent.futures.ProcessPoolExecutor(max_workers=cfg.n_procs)` で並列に回す。各ワーカーでは `OMP_NUM_THREADS=1` にする。`meta["pyxtal"]` は pickle できない可能性があるので、送る前に外す。

- [ ] **Step 1: 失敗するテストを書く**

```python
pytest.importorskip("tblite")

def test_relax_candidate_gfn2_benzene_p1_short():
    mol = load_molecule("data/benzene.toml")
    c = generate_random(mol, 2, 2, n=1, rng=np.random.default_rng(3))[0]
    out = relax_candidate(c, mol, "GFN2-xTB", fmax=0.1, steps=5, max_mol_rmsd=0.5)
    assert out.status in {"relaxed", "rejected_deformed"}
    assert out.energy_per_mol is not None and out.density > 0

def test_relax_many_survives_failure(monkeypatch):  # Review Focus 1
    mol = load_molecule("data/benzene.toml")
    cands = generate_random(mol, 2, 2, n=2, rng=np.random.default_rng(3))
    cands[0].atoms.positions[0] = cands[0].atoms.positions[1]  # coincident atoms -> calculator error
    cfg = replace(default_cfg(), n_procs=1, relax_steps=2)
    out = relax_many(cands, mol, cfg, method="GFN2-xTB")
    assert [c.status for c in out].count("failed") >= 1 and len(out) == 2

def test_make_calculator_rejects_unknown():
    with pytest.raises(ValueError):
        make_calculator("GFN-FF")
```

`default_cfg()` は Task 3 で作った `tests/csp/conftest.py` のもの。

- [ ] **Step 2: 失敗を確認する**

Run: `pytest tests/csp/test_relax.py -v`
Expected: FAIL

- [ ] **Step 3: 実装する**

- [ ] **Step 4: テストを通す**

Run: `pytest tests/csp/test_relax.py -v`
Expected: PASS（Mac で 1 分以内）

- [ ] **Step 5: コミットする**

```bash
git add source/csp/relax.py tests/csp/test_relax.py
git commit -m "feat: relax candidates with GFN2-xTB in parallel"
```

---

### Task 9: 重複除去、全体の流れ、出力、CLI

**Files:**
- Create: `source/csp/dedupe.py`、`source/csp/pipeline.py`、`source/csp/cli.py`
- Test: `tests/csp/test_pipeline.py`

**Interfaces:**
- Consumes: ここまでの全部
- Produces:
  - `deduplicate(cands: list[Candidate], ltol: float = 0.2, stol: float = 0.3, angle_tol: float = 5.0) -> list[Candidate]`。`pymatgen.analysis.structure_matcher.StructureMatcher` を使う（`pymatgen.io.ase.AseAtomsAdaptor` で変換）。エネルギーの低い方を残し、`meta["duplicates"]` に、まとめた候補の由来を足す。
  - `run_csp(cfg: CSPConfig) -> list[Candidate]`。手順は次の通り。
    1. 生成：空間群ごとに `n_candidates_per_sg` を `generators` の割合で分ける。`voronoi` は `random` で作った骨格を `orient_candidate` に通す。`seed` は `seed_structures` ごとに全件。
    2. `prescreen`
    3. `relax_many`
    4. `status == "relaxed"` のものだけを `deduplicate` に通す
    5. `energy_per_mol` の昇順に順位を付ける
    6. 出力
  - 出力は次の3つ。
    - `cfg.out_dir / "candidates.csv"`：列は `rank,id,energy_per_mol_eV,rel_energy_kJmol,density_gcm3,space_group_in,space_group_out,origin,max_mol_rmsd,status`。`space_group_out` は `spglib`（symprec=0.1）で判定する。
    - `cfg.out_dir / f"{rank:03d}_{id}.cif"`：上位 `min(50, n)` 件。
    - `cfg.out_dir / "run.log"`：各段の件数と経過時間。
  - 時間の目安：生成とふるい分けが終わった時点で、それまでの1件あたりの緩和時間の見込みと `time_budget_min` から、緩和する件数を `n_prescreen_keep` 以下に切り詰めてよい。その場合はログに書く。見込みは最初の `n_procs` 件の実測で出す。
  - `cli.main(argv: list[str] | None = None) -> int`。使い方は `molecular-modeling-csp [config]`。`config` の既定は `input/csp.toml`。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_deduplicate_merges_identical():
    c = generate_random(load_molecule("data/benzene.toml"), 2, 2, n=1, rng=np.random.default_rng(0))[0]
    a, b = replace(c, energy_per_mol=-1.0), replace(c, atoms=c.atoms.copy(), energy_per_mol=-0.5)
    assert len(deduplicate([a, b])) == 1

def test_run_csp_benzene_end_to_end(tmp_path):
    cfg = replace(default_cfg(), space_groups=[2], n_candidates_per_sg=6, n_prescreen_keep=2,
                  relax_steps=3, n_procs=1, out_dir=tmp_path)
    out = run_csp(cfg)
    rows = (tmp_path / "candidates.csv").read_text().splitlines()
    assert rows[0].startswith("rank,id,energy_per_mol_eV")
    assert len(rows) - 1 == len(out) and len(out) >= 1
    assert {c.origin for c in out} <= {"random", "voronoi"}
    assert (tmp_path / "run.log").exists()

def test_run_csp_handles_zero_survivors(tmp_path):  # Review Focus 5
    cfg = replace(default_cfg(), space_groups=[2], n_candidates_per_sg=3, packing_range=(0.99, 1.0), out_dir=tmp_path)
    assert run_csp(cfg) == []
    assert (tmp_path / "candidates.csv").read_text().splitlines() == [CSV_HEADER]

def test_cli_returns_zero(tmp_path):
    cfg_path = tmp_path / "csp.toml"
    cfg_path.write_text(
        'molecule = "data/benzene.toml"\nmode = "informed"\nspace_groups = [2]\n'
        f'n_candidates_per_sg = 4\nn_prescreen_keep = 1\nrelax_steps = 2\nn_procs = 1\nout_dir = "{tmp_path / "out"}"\n'
    )
    assert cli.main([str(cfg_path)]) == 0
    assert (tmp_path / "out" / "candidates.csv").exists()
```

`CSV_HEADER` は `pipeline.py` のモジュール定数で、テストはそれを import する。

- [ ] **Step 2: 失敗を確認する**

Run: `pytest tests/csp/test_pipeline.py -v`
Expected: FAIL

- [ ] **Step 3: 実装する**

- [ ] **Step 4: テストを通す**

Run: `pytest tests/csp/test_pipeline.py -v`
Expected: PASS

- [ ] **Step 5: コミットする**

```bash
git add source/csp/dedupe.py source/csp/pipeline.py source/csp/cli.py tests/csp/test_pipeline.py
git commit -m "feat: add CSP pipeline, outputs and CLI"
```

---

### Task 10: 答え合わせ（RMSD₁₅）

**Files:**
- Create: `source/csp/validate.py`
- Test: `tests/csp/test_validate.py`

**Interfaces:**
- Consumes: `match_and_align`、`read_seed_molecules`（実験 CIF の読み込み）、`candidates.csv` と CIF の出力
- Produces:
  - `rmsd_cluster(a: ase.Atoms, b: ase.Atoms, n_mol: int = 15) -> tuple[float, int]`。戻り値は（RMSD、一致した分子数）。手順は次の通り。
    1. a の分子0を中心に、重心が近い順に、周期像を含めて `n_mol` 個の分子を集める。
    2. b でも、各分子を中心にして同じように集める。
    3. 中心の分子同士を `match_and_align`（重原子のみ）で重ね、その変換をクラスター全体に当てはめる。
    4. クラスターの分子同士を、重心の距離で `linear_sum_assignment` を使って対応付ける。
    5. 分子内の原子の対応は、中心の分子で得た `perm` を使う。
    6. 重原子の RMSD を出し、b の中心の選び方の中で最小のものを返す。
    - 一致した分子数は、重心のずれが 2.0 Å 未満の分子の数。
  - `validate_run(out_dir: Path, experimental: list[Path], n_top: int = 20, rmsd_cut: float = 0.8, flip: tuple[str, str] | None = ("O", "S")) -> ValidationReport`。各実験構造について、次のことを調べる。
    - 全候補（上位に限らない）との RMSD₁₅
    - 上位 `n_top` 件の中に、`rmsd < rmsd_cut` かつ一致した分子数が 15 のものがあるか（あれば成功）
    - 一致した候補の順位と由来
    - `flip` を指定した場合は、実験構造の O と S を入れ替えた版とも比べ、「反転一致」として別に記録する（成功には数えない）
  - `ValidationReport` は dataclass。`to_markdown() -> str` を持つ。
  - CLI は `python -m molecular_modeling.csp.validate <out_dir> <exp.cif>...`。結果のマークダウンを標準出力に出し、`out_dir / "validation.md"` にも書く。

- [ ] **Step 1: 失敗するテストを書く**

```python
def test_rmsd_cluster_zero_for_same_structure_and_rigid_motion():
    c = generate_random(load_molecule("data/benzene.toml"), 14, 4, n=1, rng=np.random.default_rng(5))[0].atoms
    moved = c.copy(); moved.translate([0.3, -0.2, 0.1]); moved.wrap()
    rmsd, n = rmsd_cluster(c, moved)
    assert rmsd < 1e-3 and n == 15

def test_rmsd_cluster_large_for_different_packing():
    mol = load_molecule("data/benzene.toml")
    a = generate_random(mol, 14, 4, n=1, rng=np.random.default_rng(5))[0].atoms
    b = generate_random(mol, 2, 2, n=1, rng=np.random.default_rng(6))[0].atoms
    rmsd, n = rmsd_cluster(a, b)
    assert rmsd > 0.8 or n < 15

def test_validate_run_reports_rank_and_origin(tmp_path):
    cfg = replace(default_cfg(), space_groups=[2], n_candidates_per_sg=6, n_prescreen_keep=2,
                  relax_steps=3, n_procs=1, out_dir=tmp_path)
    out = run_csp(cfg)
    exp = tmp_path / "exp.cif"
    shutil.copy(next(tmp_path.glob("001_*.cif")), exp)
    report = validate_run(tmp_path, [exp], flip=None)
    entry = report.entries[0]
    assert report.success and entry.best_rank == 1 and entry.best_origin == out[0].origin
```

`ValidationReport` は `success: bool` と `entries: list[ValidationEntry]` を持つ。`ValidationEntry` は `experimental: str`、`best_rank: int | None`、`best_rmsd: float | None`、`best_origin: str | None`、`flip_match_rank: int | None` を持つ。

- [ ] **Step 2: 失敗を確認する**

Run: `pytest tests/csp/test_validate.py -v`
Expected: FAIL

- [ ] **Step 3: 実装する**

- [ ] **Step 4: テストを通す**

Run: `pytest tests/csp/test_validate.py -v`
Expected: PASS

- [ ] **Step 5: コミットする**

```bash
git add source/csp/validate.py tests/csp/test_validate.py
git commit -m "feat: compare CSP candidates with experimental structures"
```

---

### Task 11: 実測による既定値の調整とドキュメント

**Files:**
- Modify: `source/csp/config.py`（既定値）、`input/csp.toml`
- Create: `docs/ja/csp_guide.rst`、`docs/en/csp_guide.rst`（`docs/*/index.rst` の toctree に足す）
- Test: `tests/csp/test_pipeline.py` に `@pytest.mark.slow` のテスト `test_run_csp_time_budget_respected` を足す

**Interfaces:**
- Consumes: `run_csp`

- [ ] **Step 1: 失敗するテストを書く**

```python
@pytest.mark.slow
def test_run_csp_time_budget_respected(tmp_path):
    cfg = replace(default_cfg(), space_groups=[2, 14], n_candidates_per_sg=100, time_budget_min=3, out_dir=tmp_path)
    t0 = time.monotonic(); run_csp(cfg)
    assert time.monotonic() - t0 < 3 * 60 * 1.2
```

- [ ] **Step 2: 実行して結果を見る**

Run: `pytest tests/csp/test_pipeline.py::test_run_csp_time_budget_respected -v -m slow`
Expected: Task 9 の時間の目安のロジックが正しければ PASS。FAIL なら、Task 9 の切り詰めを直す。

- [ ] **Step 3: OS-2SMe 程度の大きさで実測する**

加藤さんの OS-2SMe 分子がまだ届いていなければ、重原子が 15〜20 個程度の代わりの分子（たとえば `data/naphthalene.toml`）で測る。次の3つを測る。
- GFN2-xTB の緩和1件あたりの時間（`relax_steps` = 150、Z = 4）
- PyXtal の生成1件あたりの時間
- 剛体 LJ の1件あたりの時間

1時間・Mac の全コアで収まるように、`n_candidates_per_sg`、`n_prescreen_keep`、`relax_steps` の既定値を決めて `config.py` に書く。測った値は `run.log` の書式でガイドに載せる。

- [ ] **Step 4: ガイドを書く（日英）**

ガイドに書くこと:
- インストール（`pip install -e ".[csp]"`）
- `input/csp.toml` の各キー
- 2つのモード（blind、informed）
- 出力の読み方
- 答え合わせのやり方
- 既存の Voronoi パイプラインとの関係（Voronoi は「向きの初期値」として使うようになったこと）

Run: `SPHINX_LANG=ja make -C docs html && SPHINX_LANG=en make -C docs html`
Expected: 警告なしでビルドが通る

- [ ] **Step 5: 全体の確認とコミット**

Run: `pytest -m "not slow" && ruff check source/ tests/ && mypy source/ --ignore-missing-imports`
Expected: すべて PASS、エラーなし

```bash
git add source/csp/config.py input/csp.toml docs/ja/csp_guide.rst docs/en/csp_guide.rst docs/ja/index.rst docs/en/index.rst tests/csp/test_pipeline.py
git commit -m "docs: add CSP guide and tune defaults from measured timings"
```

---

## 計画の範囲外（設計書の通り）

- OS-2SMe での本番の予測と答え合わせ。加藤さんの分子構造と実験 CIF が届いてから、この計画の CLI で行う。
- flip disorder の再現、DFT、スパコン。
