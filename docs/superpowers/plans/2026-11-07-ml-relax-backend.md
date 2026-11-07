# CSP 緩和の計算方法選択（OrbMol-v2）実装計画

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** CSP の緩和で `relax_method` により GFN2-xTB（既定）／GFN1-xTB／OrbMol-v2 を選べるようにし、既存の時間予算・失敗処理・出力の流れをそのまま保つ。

**Architecture:** 計算器の生成・動作確認・バージョン情報を新モジュール `source/csp/backends.py` にまとめる（親プロセスでは torch を import しない）。`relax.py` に、spawn で起動し準備ジョブでモデルを読み込ませる再利用可能なプール `WorkerPool` を足し、`_relax_with_budget` が 1 回の緩和段階を通してそれを使い回す。設定・CLI・`run.log`・`candidates.csv` に方法名と出どころを通す。

**Tech Stack:** Python 3.11+, ASE, tblite, orb-models 0.7.x（PyTorch, NVIDIA Warp）, concurrent.futures / multiprocessing（spawn）, pytest。

**Spec:** `docs/superpowers/specs/2026-11-07-ml-relax-backend-design.md`（利用者承認済み 2026-11-07。既定値は `"GFN2-xTB"` のまま）

## Global Constraints

- 入力ファイルで選べる値は `"GFN2-xTB"`, `"GFN1-xTB"`, `"OrbMol-v2"` の3つだけ。既定 `"GFN2-xTB"`。`EMT` 系（テスト用）は `run_csp(..., method=...)` からのみ。
- extra: `csp-ml = ["molecular_modeling[csp]", "orb-models>=0.7,<0.8"]`。足りないときの案内は `pip install -e ".[csp-ml]"`。
- 親プロセスでは torch・warp・orb_models を import しない（`importlib.util.find_spec` で存在確認のみ）。
- すべてのワーカープール（並列・単独・動作確認）は `multiprocessing.get_context("spawn")` で作る。
- OrbMol-v2: `orbmol_v2(device="cpu")` と `ORBCalculator(model, adapter, device="cpu")`。緩和前に `atoms.info["charge"] = 0`, `atoms.info["spin"] = 1`。
- ワーカー初期化で `OMP_NUM_THREADS=1`、torch が読み込まれていれば `torch.set_num_threads(1)`。
- 動作確認のタイムアウトは 600 秒。失敗は `RelaxBackendError(RuntimeError)`。`cli.EXPECTED_ERRORS` に足すのはこの型だけ。
- 動作確認は OrbMol-v2 のときだけ、出力ディレクトリに触る前に行う。
- 出どころの値は `default` / `config` / `cli` / `argument` の4つ。
- `candidates.csv` の末尾に `relax_method` 列。値は実際に使った方法名。
- `relax_steps`・`relax_fmax`・事前ふるいは方法によらず共通。
- 既定（xTB）の利用者の結果を変えない：xTB と EMT の計算器はワーカー内でキャッシュせず、候補ごとに新しく作る（今と同じ）。
- 重いモデルを読むテストは `@pytest.mark.slow` かつ `orb_models` が無ければスキップ。CI の通常テストの時間を増やさない。
- コミットメッセージ・PR 本文に社内ツール名を書かない。末尾に `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>` と `Claude-Session: https://claude.ai/code/session_01AfapaMtNCZsX6hYifZCkkx`。

## Review Focus

1. **OrbMol-v2 を選んだのに重みが取れない（オフライン）** → 出力ディレクトリに触る前に `error:` 1 行（方法名・原因・キャッシュ場所・対処）で止まる。Task 4 の `test_run_csp_preflight_failure_stops_before_out_dir` と `test_cli_reports_backend_error_in_one_line`。
2. **`n_procs = 1` で何ラウンドも回す** → モデル読み込みはプール 1 個分だけ（ラウンドごとに作り直さない）。Task 3 の `test_budget_reuses_one_pool_across_rounds`（`n_procs=1` を含む）。
3. **途中のラウンドでワーカーが落ちる** → 落とした候補だけ `failed`、次のラウンドの前にプールを作り直して準備し直し、以後のラウンドは普通に進む。Task 3 の `test_budget_recovers_after_worker_crash_mid_run`。
4. **準備（モデル読み込み）が長い** → 経過時間には入るが 1 ラウンドの見積もりには入らない。予算が準備だけで尽きても最初のラウンドは走る。Task 3 の `test_budget_counts_warmup_as_elapsed_not_round_time` と `test_budget_first_round_runs_even_if_warmup_eats_budget`。
5. **`--relax-method` と設定ファイルの値が食い違う** → CLI が勝ち、`run.log` に `cli` と書かれ、CSV の列も CLI の値。Task 2 の `test_cli_relax_method_overrides_config`。

---

## File Structure

| ファイル | 役割 |
|---|---|
| Create `source/csp/backends.py` | 方法名の定数、`make_calculator`、`get_calculator`（ワーカー内キャッシュ）、`RelaxBackendError`、`preflight`、`backend_info`。トップレベルで torch を import しない |
| Modify `source/csp/relax.py` | `make_calculator`/`METHODS` を backends から再公開。spawn 化、`_init_worker(method)`、`WorkerPool`、`relax_many(..., pool=None)`、charge/spin |
| Modify `source/csp/deps.py` | `require(..., hint=INSTALL_HINT)` |
| Modify `source/csp/config.py` | `relax_method`、`relax_method_source`、`CONFIG_RELAX_METHODS` と検査 |
| Modify `source/csp/pipeline.py` | `check_dependencies` の方法別化、`run_csp(cfg, method=None)`、動作確認、`run.log` の記録、`WorkerPool` を使う `_relax_with_budget`、CSV 列 |
| Modify `source/csp/cli.py` | `--relax-method`、`EXPECTED_ERRORS` に `RelaxBackendError` |
| Modify `pyproject.toml`, `uv.lock` | `csp-ml` extra |
| Tests | `tests/csp/test_backends.py`（新規）、`test_relax.py`、`test_config.py`、`test_pipeline.py` |
| Docs | `docs/{en,ja}/csp_guide.rst`, `installation.rst`, `cli_reference.rst`, `README.md`, `README_ja.md`, `input/csp.toml`, `CLAUDE.md` |

---

### Task 1: 計算器モジュール `backends.py` と依存チェック

**Files:**
- Create: `source/csp/backends.py`, `tests/csp/test_backends.py`
- Modify: `source/csp/relax.py:22-52`（`METHODS`・`make_calculator` を移して再公開）, `source/csp/deps.py`, `source/csp/pipeline.py:42-57`, `pyproject.toml`, `uv.lock`

**Interfaces:**
- Produces（`molecular_modeling.csp.backends`）:
  - `CONFIG_METHODS: tuple[str, ...]` は `config.CONFIG_RELAX_METHODS` と同じ値 `("GFN2-xTB", "GFN1-xTB", "OrbMol-v2")`（定義は config.py 側に置き、backends が import する。Task 2 で config に追加するので、Task 1 では config.py に定数だけ先に置く）
  - `TEST_METHODS = ("EMT", "EMT-CRASH-ON-AU", "EMT-BROKEN", "EMT-SLOW")`
  - `METHODS = CONFIG_METHODS + TEST_METHODS`
  - `ML_METHODS = frozenset({"OrbMol-v2"})`（torch を使い、ワーカー内キャッシュと動作確認の対象）
  - `make_calculator(method: str) -> Calculator`
  - `class RelaxBackendError(RuntimeError)`
  - `ML_INSTALL_HINT = 'install the machine-learned relaxation dependencies with: pip install -e ".[csp-ml]"'`
- Produces（`deps.require`）: `require(modules, attributes=None, purpose="this command", hint: str = INSTALL_HINT) -> None`
- Produces（`pipeline.check_dependencies(method: str = "GFN2-xTB") -> None`）: 既存の挙動 + OrbMol-v2 の分岐

- [ ] **Step 1: 失敗するテストを書く**（`tests/csp/test_backends.py`）

```python
def test_methods_lists():
    assert CONFIG_METHODS == ("GFN2-xTB", "GFN1-xTB", "OrbMol-v2")
    assert set(TEST_METHODS) <= set(METHODS) and "OrbMol-v2" in ML_METHODS

def test_make_calculator_test_methods():
    assert make_calculator("EMT").__class__.__name__ == "EMT"
    with pytest.raises(RuntimeError, match="test backend is broken"):
        make_calculator("EMT-BROKEN")
    with pytest.raises(ValueError, match="GFN-FF"):
        make_calculator("GFN-FF")

def test_importing_backends_does_not_import_torch():
    # 新しいインタープリタで確かめる（他テストが torch を読んでいても影響しない）
    code = "import sys, molecular_modeling.csp.backends, molecular_modeling.csp.pipeline; print('torch' in sys.modules)"
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False"

def test_check_dependencies_orbmol_names_ml_extra(monkeypatch):
    real = importlib.util.find_spec
    monkeypatch.setattr(importlib.util, "find_spec", lambda n, *a: None if n in ("orb_models", "torch") else real(n, *a))
    with pytest.raises(ImportError, match=r"orb_models.*torch.*\.\[csp-ml\]"):
        check_dependencies("OrbMol-v2")

def test_check_dependencies_orbmol_does_not_need_tblite(monkeypatch):
    # tblite.ase の import を失敗させても、orb_models/torch が「ある」なら通る
    ...  # fake import_module で tblite* を ImportError、fake find_spec で orb_models/torch を非 None
    check_dependencies("OrbMol-v2")
```

`EMT-SLOW` は `make_calculator` 内で `time.sleep(30)` してから EMT を返す（Task 4 のタイムアウト試験用）。`EMT-BROKEN` は `RuntimeError("test backend is broken")` を投げる。

- [ ] **Step 2: 失敗を確認**

Run: `pytest tests/csp/test_backends.py -v`
Expected: FAIL（`ModuleNotFoundError: molecular_modeling.csp.backends`）

- [ ] **Step 3: 実装**
  - `config.py` に `CONFIG_RELAX_METHODS = ("GFN2-xTB", "GFN1-xTB", "OrbMol-v2")` を置く（検査は Task 2）。
  - `backends.py` に上の定数・`RelaxBackendError`・`make_calculator` を置く。OrbMol-v2 分岐は関数内で `from orb_models.forcefield import pretrained` と `from orb_models.forcefield.inference.calculator import ORBCalculator` を import。
  - `relax.py` は `from .backends import METHODS, make_calculator` で再公開（既存テストの import を保つ）。旧 `make_calculator` の本体は削除。
  - `deps.require` に `hint` 引数を足し、メッセージ末尾を `hint` にする。
  - `check_dependencies`: `OrbMol-v2` は `require(REQUIRED_PACKAGES, purpose=..., hint=ML_INSTALL_HINT)` に加え、`find_spec` が `None` の `orb_models`・`torch` を集めて、同じ書式の `ImportError`（名前を列挙し `ML_INSTALL_HINT` で終わる）を投げる。GFN 系と EMT 系は今のまま。
  - `pyproject.toml` に `csp-ml` extra を足し、`uv lock` で `uv.lock` を更新する。

- [ ] **Step 4: 通ることを確認**

Run: `pytest tests/csp -m "not slow" -q`
Expected: 全 PASS（既存の `test_make_calculator_rejects_unknown`, `test_unknown_method_rejected_before_anything_runs` を含む）

- [ ] **Step 5: lint とコミット**

Run: `ruff check source/ tests/ && ruff format --check source/ tests/ && mypy source/ --ignore-missing-imports`

```bash
git add source/csp/backends.py source/csp/relax.py source/csp/deps.py source/csp/pipeline.py source/csp/config.py tests/csp/test_backends.py pyproject.toml uv.lock
git commit -m "Add a calculator module with an optional OrbMol-v2 backend"
```

---

### Task 2: 設定・CLI・CSV 列・出どころの記録

**Files:**
- Modify: `source/csp/config.py`, `source/csp/cli.py`, `source/csp/pipeline.py`（`run_csp`, `_run_stages`, `_write_outputs`, `CSV_HEADER`）
- Test: `tests/csp/test_config.py`, `tests/csp/test_pipeline.py`

**Interfaces:**
- Consumes: `CONFIG_RELAX_METHODS`（Task 1）
- Produces:
  - `CSPConfig.relax_method: str = "GFN2-xTB"`
  - `CSPConfig.relax_method_source: str = field(default="default", repr=False, compare=False)`（`_NOT_INPUT_KEYS` に追加）
  - `run_csp(cfg: CSPConfig, method: str | None = None) -> list[Candidate]`
  - `pipeline.resolve_method(cfg: CSPConfig, method: str | None) -> tuple[str, str]`（方法名, 出どころ）
  - `CSV_HEADER` の末尾に `,relax_method`
  - `_write_outputs(cfg, ranked, method: str)`
  - `run.log` の行: `relax method %s (from %s)`（Task 4 でバージョン等を同じ行の後ろに足す）

- [ ] **Step 1: 失敗するテストを書く**

`test_config.py`:
```python
def test_relax_method_default_and_source(tmp_path):
    cfg = load_csp_config(_toml(tmp_path, 'molecule = "data/benzene.toml"\n'))
    assert (cfg.relax_method, cfg.relax_method_source) == ("GFN2-xTB", "default")

def test_relax_method_from_file(tmp_path):
    cfg = load_csp_config(_toml(tmp_path, 'molecule = "data/benzene.toml"\nrelax_method = "OrbMol-v2"\n'))
    assert (cfg.relax_method, cfg.relax_method_source) == ("OrbMol-v2", "config")

@pytest.mark.parametrize("bad", ['"EMT"', '"orbmol-v2"', "1"])
def test_relax_method_rejects_other_values(tmp_path, bad):
    with pytest.raises(ValueError, match="relax_method must be one of 'GFN2-xTB', 'GFN1-xTB', 'OrbMol-v2'"):
        load_csp_config(_toml(tmp_path, f'molecule = "data/benzene.toml"\nrelax_method = {bad}\n'))

def test_relax_method_source_is_not_an_input_key(tmp_path):
    with pytest.raises(ValueError, match="unknown key"):
        load_csp_config(_toml(tmp_path, 'molecule = "data/benzene.toml"\nrelax_method_source = "cli"\n'))
```
（`_toml` は test_config.py の既存ヘルパーに合わせる。無ければファイルを書いて Path を返す小関数を足す。）

`test_pipeline.py`:
```python
def test_csv_has_relax_method_column(tmp_path, monkeypatch):
    # EMT で 1 件以上残る小さな実行（既存の end_to_end と同じ設定）を method="EMT" で
    ...
    header, first = (tmp_path / "candidates.csv").read_text().splitlines()[:2]
    assert header.endswith(",relax_method") and first.endswith(",EMT")

@pytest.mark.parametrize(("method", "source"), [(None, "default"), ("EMT", "argument")])
def test_run_log_records_method_and_source(tmp_path, monkeypatch, method, source):
    monkeypatch.setattr(pipeline, "_generate", lambda *a, **k: [])
    monkeypatch.setattr(pipeline, "check_dependencies", lambda m: None)  # xTB 不要
    run_csp(replace(default_cfg(), out_dir=tmp_path), method=method)
    expect = "GFN2-xTB" if method is None else "EMT"
    assert f"relax method {expect} (from {source})" in (tmp_path / "run.log").read_text()

def test_cli_relax_method_overrides_config(tmp_path, monkeypatch):
    seen = {}
    monkeypatch.setattr(cli, "run_csp", lambda cfg: seen.update(m=cfg.relax_method, s=cfg.relax_method_source))
    p = tmp_path / "csp.toml"
    p.write_text('molecule = "data/benzene.toml"\nrelax_method = "GFN2-xTB"\n')
    assert cli.main([str(p), "--relax-method", "GFN1-xTB"]) == 0
    assert seen == {"m": "GFN1-xTB", "s": "cli"}

def test_cli_rejects_unknown_relax_method(tmp_path, capsys):
    with pytest.raises(SystemExit):
        cli.main(["x.toml", "--relax-method", "EMT"])
```

- [ ] **Step 2: 失敗を確認**

Run: `pytest tests/csp/test_config.py tests/csp/test_pipeline.py -k "relax_method or csv_has" -v`
Expected: FAIL

- [ ] **Step 3: 実装**
  - ローダー: `relax_method` が文字列で `CONFIG_RELAX_METHODS` に含まれなければ `ValueError("relax_method must be one of 'GFN2-xTB', 'GFN1-xTB', 'OrbMol-v2', got <repr>")`。ファイルにあれば `relax_method_source = "config"`。
  - CLI: `--relax-method` を `choices=CONFIG_RELAX_METHODS` で足し、指定時は `cfg.relax_method` を書き換えて `relax_method_source = "cli"`。`run_csp(cfg)` を呼ぶ（引数で上書きしない）。
  - `resolve_method`: `method` が `None` なら `(cfg.relax_method, cfg.relax_method_source)`、そうでなければ `(method, "argument")`。`run_csp` の先頭で使い、以降はこの方法名だけを使う。`cfg` は書き換えない。
  - `run.log` を開いた直後（`config_warnings` を書く場所の直後）に `logger.info("relax method %s (from %s)", method, source)`。
  - `CSV_HEADER` 末尾に `relax_method`、各行の末尾に方法名。
  - 既存テストの `fake_relax_many`／`_old_results` などで新しい列や引数に合わせて直す必要があれば直す。

- [ ] **Step 4: 通ることを確認**

Run: `pytest tests/csp -m "not slow" -q`
Expected: 全 PASS（`test_validate.py` を含む。validate は `rank,id,origin` しか読まない）

- [ ] **Step 5: lint とコミット**

```bash
git commit -am "Let CSP runs choose the relaxation method in the config or on the command line"
```

---

### Task 3: spawn・再利用プール・準備ジョブ・時間予算

**Files:**
- Modify: `source/csp/backends.py`（`get_calculator`）, `source/csp/relax.py`（`_init_worker`, `_worker`, `_run_parallel`, `_run_isolated`, `relax_candidate`, `relax_many`, 新 `WorkerPool`）, `source/csp/pipeline.py:183-213`
- Test: `tests/csp/test_relax.py`, `tests/csp/test_pipeline.py`

**Interfaces:**
- Consumes: `make_calculator`, `ML_METHODS`, `RelaxBackendError`（Task 1）
- Produces:
  - `backends.get_calculator(method: str) -> Calculator`: `ML_METHODS` の方法はモジュール変数 `_CACHE: dict[str, Calculator]` に保持し、返す前に `calc.reset()`。それ以外は毎回 `make_calculator(method)`。
  - `relax._init_worker(method: str) -> None`: `OMP_NUM_THREADS=1` → `get_calculator(method)` → `"torch" in sys.modules` なら `torch.set_num_threads(1)`。
  - `class WorkerPool`（`relax.py`）:
    - `__init__(self, workers: int, method: str)`（まだ起動しない）
    - `ensure(self) -> float`: プールが無いか壊れていれば、`ProcessPoolExecutor(max_workers=workers, mp_context=get_context("spawn"), initializer=_init_worker, initargs=(method,))` を作り、`_warmup` ジョブを `workers` 個投げて全部待ち、かかった秒数を返して `logger.info("worker pool: %d workers ready in %.1f s (%s)", ...)`。健全なら `0.0`。準備で `BrokenProcessPool` などが出たら閉じて `RelaxBackendError("relaxation workers could not start (<method>): <type>: <msg>")`。
    - `executor` プロパティ、`mark_broken(self) -> None`（shutdown して捨てる）、`close(self) -> None`、`__enter__`/`__exit__`。
  - `relax_many(cands, ref_mol, cfg, method="GFN2-xTB", pool: WorkerPool | None = None) -> list[Candidate]`: `pool` があれば `n_procs` によらず並列段は `pool.executor` で回し（`with` で閉じない）、壊れたら `pool.mark_broken()`。`pool` が無ければ今と同じ（自前のプールを作って閉じる）。単独再実行（`_run_isolated`）は今の仕組みのまま、ただし spawn・`initargs=(method,)`。
  - `_relax_with_budget(cands, mol, cfg, method, t_start)` は署名を変えず、内部で `with WorkerPool(max(cfg.n_procs, 1), method) as pool:` を作り、各ラウンドの前に `pool.ensure()` を呼んで（ラウンド計時 `t0` の**前**）、`relax_many(batch, mol, cfg, method, pool=pool)`。
  - `relax_candidate` は `get_calculator(method)` を使い、`method in ML_METHODS` なら `atoms.info.update(charge=0, spin=1)`。

- [ ] **Step 1: 失敗するテストを書く**

`test_relax.py`（`_ControlledPool.__init__` を `(self, max_workers=None, **kwargs)` にして `kwargs` を記録する形に直す）:
```python
def test_all_pools_use_spawn(monkeypatch):
    # _ControlledPool.created に kwargs を積み、crash={1} で並列・単独の両方を通す
    ...
    assert all(k["mp_context"].get_start_method() == "spawn" for k in _ControlledPool.created)
    assert all(k["initargs"] == ("EMT",) for k in _ControlledPool.created)

def test_relax_many_uses_given_pool_and_leaves_it_open():
    mol = Atoms("Cu2", positions=[[0, 0, 0], [0, 0, 2.5]])
    cfg = replace(default_cfg(), n_procs=2, relax_steps=2, max_mol_rmsd=5.0)
    with WorkerPool(2, "EMT") as pool:
        pool.ensure()
        ex = pool.executor
        relax_many(_emt_cands(2), mol, cfg, "EMT", pool=pool)
        relax_many(_emt_cands(2), mol, cfg, "EMT", pool=pool)
        assert pool.executor is ex and pool.ensure() == 0.0

def test_worker_pool_rebuilds_after_crash():
    ...  # EMT-CRASH-ON-AU で gold={0} の 2 件 → pool が壊れ、次の ensure() が > 0 を返し新しい executor

def test_worker_pool_start_failure_is_backend_error():
    with WorkerPool(1, "EMT-BROKEN") as pool, pytest.raises(RelaxBackendError, match="EMT-BROKEN"):
        pool.ensure()

def test_get_calculator_caches_only_ml_methods(monkeypatch):
    monkeypatch.setattr(backends, "ML_METHODS", frozenset({"EMT"}))
    assert backends.get_calculator("EMT") is backends.get_calculator("EMT")
    monkeypatch.setattr(backends, "ML_METHODS", frozenset())
    backends._CACHE.clear()
    assert backends.get_calculator("EMT") is not backends.get_calculator("EMT")
```

`test_pipeline.py` の `_run_budget` を、`pipeline.WorkerPool` を偽物（`ensure()` が初回と `mark_broken` 後に `clock.now += warmup_s` して `warmup_s` を返す、生成回数を数える）に差し替える形に直し、`fake_relax_many(cands, mol, cfg, method, pool=None)` にする。既存 4 本の期待値はそのまま通ること（`warmup_s=0` 既定）。追加:
```python
@pytest.mark.parametrize("n_procs", [1, 4])
def test_budget_reuses_one_pool_across_rounds(monkeypatch, n_procs):
    _, sizes, _, pools = _run_budget(monkeypatch, budget_min=60, n_cands=10, n_procs=n_procs)
    assert len(sizes) > 1 and pools.created == 1

def test_budget_counts_warmup_as_elapsed_not_round_time(monkeypatch):
    # 予算 540 s、準備 200 s、1 ラウンド 100 s → 200 + 3*100 = 500 で残り 40 s < 100 s、3 ラウンドで止まる
    done, sizes, t, _ = _run_budget(monkeypatch, budget_min=10, warmup_s=200.0)
    assert sizes == [4] * 3 and t == 500

def test_budget_first_round_runs_even_if_warmup_eats_budget(monkeypatch):
    done, sizes, _, _ = _run_budget(monkeypatch, budget_min=1, warmup_s=1000.0)
    assert sizes == [4]

def test_budget_recovers_after_worker_crash_mid_run(tmp_path):
    # 実プロセス・EMT-CRASH-ON-AU、n_procs=2、6 候補（2 ラウンド目に金を 1 つ）、十分な予算
    out = pipeline._relax_with_budget(cands, mol, cfg, "EMT-CRASH-ON-AU", time.perf_counter())
    assert [c.status for c in out].count("failed") == 1 and len(out) == 6
```
（`_run_budget` の戻り値に偽プールを足すので、既存 4 本のアンパックも直す。）

- [ ] **Step 2: 失敗を確認**

Run: `pytest tests/csp/test_relax.py tests/csp/test_pipeline.py -m "not slow" -v`
Expected: 新規テストが FAIL

- [ ] **Step 3: 実装**（上の Interfaces どおり。`_warmup()` はモジュール関数で `os.getpid()` を返すだけ。`relax.py` に `logger = logging.getLogger(__name__)` を足す）

- [ ] **Step 4: 通ることを確認**

Run: `pytest tests/csp -m "not slow" -q`
Expected: 全 PASS。所要時間が以前と大きく変わらないこと（spawn で数秒増える程度）。

- [ ] **Step 5: lint とコミット**

```bash
git commit -am "Start relaxation workers with spawn and reuse one warmed-up pool per relaxation stage"
```

---

### Task 4: 動作確認・バージョン記録・1 行エラー

**Files:**
- Modify: `source/csp/backends.py`（`preflight`, `backend_info`）, `source/csp/pipeline.py`（`run_csp`）, `source/csp/cli.py`
- Test: `tests/csp/test_backends.py`, `tests/csp/test_pipeline.py`

**Interfaces:**
- Consumes: `RelaxBackendError`, `make_calculator`, `ML_METHODS`（Task 1）, `resolve_method`（Task 2）
- Produces:
  - `PREFLIGHT_METHODS = frozenset({"OrbMol-v2"})`（`pipeline` から参照。テストで差し替える）
  - `preflight(method: str, timeout_s: float = 600.0) -> dict[str, str]`: `get_context("spawn").Pool(1)` で `_preflight_child(method)` を `apply_async(...).get(timeout_s)`、`finally` で `terminate()`/`join()`。子は計算器を作り、水 1 分子（一辺 10 Å の周期セル、`pbc=True`、`info` に charge 0/spin 1）でエネルギー・力・応力を計算し、`{"weights": ..., "cache": ...}` を返す（EMT 系は空 dict）。子の中の例外は子で捕まえて `("error", "<Type>: <msg>")` として返す。失敗・タイムアウトは `RelaxBackendError`。
  - メッセージ（1 行）: `relaxation method <method> is not usable: <Type>: <msg>; weights are cached in <cache>; on a computer without network access, run python -c "from orb_models.forcefield import pretrained; pretrained.orbmol_v2(device='cpu')" on one with access and copy the <cache> directory`。タイムアウト時の原因部分は `TimeoutError: no answer within <timeout_s> s`。`<cache>` は `os.environ.get("CACHED_PATH_CACHE_ROOT", "~/.cache/cached_path")`（親で求め、torch を import しない）。
  - `backend_info(method: str) -> str`: `importlib.metadata.version` で `tblite <ver>` / `orb-models <ver>` / `test calculator`。
  - `run.log` の行（Task 2 の行を拡張）: `relax method %s (from %s); %s` に続けて、動作確認の結果があれば `; weights %s in %s`。
  - `cli.EXPECTED_ERRORS` に `RelaxBackendError` を追加。

- [ ] **Step 1: 失敗するテストを書く**

```python
# test_backends.py
def test_preflight_ok_for_working_backend():
    assert preflight("EMT") == {}

def test_preflight_failure_is_one_line_with_advice():
    with pytest.raises(RelaxBackendError) as ei:
        preflight("EMT-BROKEN")
    msg = str(ei.value)
    assert "\n" not in msg and "EMT-BROKEN" in msg and "test backend is broken" in msg
    assert "cached_path" in msg and "pretrained.orbmol_v2" in msg

def test_preflight_times_out():
    t0 = time.monotonic()
    with pytest.raises(RelaxBackendError, match="TimeoutError"):
        preflight("EMT-SLOW", timeout_s=2)
    assert time.monotonic() - t0 < 15

# test_pipeline.py
def test_run_csp_preflight_failure_stops_before_out_dir(tmp_path, monkeypatch):
    monkeypatch.setattr(pipeline, "PREFLIGHT_METHODS", frozenset({"EMT-BROKEN"}))
    out = tmp_path / "out"
    with pytest.raises(RelaxBackendError):
        run_csp(replace(default_cfg(), out_dir=out), method="EMT-BROKEN")
    assert not out.exists()

def test_xtb_runs_skip_preflight(tmp_path, monkeypatch):
    monkeypatch.setattr(pipeline, "preflight", lambda *a, **k: pytest.fail("no preflight for xTB"))
    monkeypatch.setattr(pipeline, "_generate", lambda *a, **k: [])
    run_csp(replace(default_cfg(), out_dir=tmp_path))
    assert "tblite" in (tmp_path / "run.log").read_text()

def test_cli_reports_backend_error_in_one_line(tmp_path, monkeypatch, capsys):
    def boom(cfg):
        raise RelaxBackendError("relaxation method OrbMol-v2 is not usable: OSError: offline; ...")
    monkeypatch.setattr(cli, "run_csp", boom)
    p = tmp_path / "csp.toml"
    p.write_text('molecule = "data/benzene.toml"\n')
    assert cli.main([str(p)]) == 2
    err = capsys.readouterr().err
    assert err.startswith("error: relaxation method OrbMol-v2") and "Traceback" not in err
```

- [ ] **Step 2: 失敗を確認** — Run: `pytest tests/csp/test_backends.py tests/csp/test_pipeline.py -k "preflight or backend_error" -v` / Expected: FAIL

- [ ] **Step 3: 実装** — `run_csp` の順序: `check_single_process` → `resolve_method` → `check_dependencies(method)` → `method in PREFLIGHT_METHODS` なら `info = preflight(method)` → 分子読み込み以降（既存）→ ログを開いた直後に方法・出どころ・`backend_info`・`info` を 1 行で書く。

- [ ] **Step 4: 通ることを確認** — Run: `pytest tests/csp -m "not slow" -q` / Expected: 全 PASS

- [ ] **Step 5: lint とコミット**

```bash
git commit -am "Check the machine-learned backend before a run and report failures in one line"
```

---

### Task 5: 実モデルでの確認とドキュメント

**Files:**
- Test: `tests/csp/test_backends.py`（slow）
- Modify: `docs/{en,ja}/csp_guide.rst`, `docs/{en,ja}/installation.rst`, `docs/{en,ja}/cli_reference.rst`, `README.md`, `README_ja.md`, `input/csp.toml`, `CLAUDE.md`

- [ ] **Step 1: `.venv` に `csp-ml` を入れる** — Run: `pip install -e ".[csp-ml]"` / Expected: `orb-models 0.7.x` と torch が入る。初回の重み取得にはネットが要る。

- [ ] **Step 2: slow テストを書く**

```python
needs_orb = pytest.mark.skipif(importlib.util.find_spec("orb_models") is None, reason="orb-models not installed")

@pytest.mark.slow
@needs_orb
def test_preflight_orbmol_real():
    info = preflight("OrbMol-v2")
    assert info["weights"] and info["cache"]

@pytest.mark.slow
@needs_orb
def test_run_csp_orbmol_small(tmp_path):
    cfg = replace(default_cfg(), n_candidates_per_sg=4, n_prescreen_keep=2, relax_steps=3, n_procs=2, out_dir=tmp_path)
    out = run_csp(cfg, method="OrbMol-v2")
    log = (tmp_path / "run.log").read_text()
    assert "relax method OrbMol-v2 (from argument); orb-models 0.7" in log and "worker pool" in log
    assert all(r.endswith(",OrbMol-v2") for r in (tmp_path / "candidates.csv").read_text().splitlines()[1:])
    assert len(out) >= 1
```

- [ ] **Step 3: 実行** — Run: `pytest tests/csp -m slow -k orbmol -v` / Expected: PASS。あわせて `pytest tests/csp -q`（slow を含む全体）も PASS。

- [ ] **Step 4: ドキュメント**（日英で同じ内容）
  - `csp_guide.rst`: 新節「計算方法の選び方」（`relax_method` の 3 値、OrbMol-v2 の入れ方 `pip install -e ".[csp-ml]"`、初回の重み取得とキャッシュ場所、オフラインでの事前取得の 1 行コマンドとキャッシュをディレクトリごと持ってくる方法、Se を含む分子などでの OrbMol-v2 の勧め（2026-11-06 の比較試験に基づく、ただし少数系での試験である旨）、動作確認と 1 行エラー、`run.log` の記録）。「GFN2-xTB」の決め打ち表現（概要、処理の流れ、`n_prescreen_keep` の説明、`energy_per_mol_eV`、限界）を方法名に依らない表現にし、「エネルギーの絶対値は方法ごとに基準が違い、比べてよいのは同じ方法の中だけ」の注意と `relax_method` 列を CSV の表に足す。設定表に `relax_method` 行。
  - `installation.rst`: `csp-ml` extra（約 1 GB、Apache 2.0 の依存）。
  - `cli_reference.rst`: `--relax-method`（設定ファイルより優先）と、`error:` になる場合に動作確認の失敗を足す。
  - `README.md`/`README_ja.md`: CSP の行を「GFN2-xTB（既定）または OrbMol-v2 で緩和」に。
  - `input/csp.toml`: `relax_method = "GFN2-xTB"` の行とコメント（OrbMol-v2 の勧め・extra）、`relax_steps` のコメントから「GFN2-xTB」を外す。
  - `CLAUDE.md`: CSP の段落の `relax (GFN2-xTB via tblite, ...)` を方法選択に合わせ、`backends.py` を挙げる。

- [ ] **Step 5: ドキュメントのビルド確認** — Run: `SPHINX_LANG=en make -C docs html && SPHINX_LANG=ja make -C docs html` / Expected: 新しい警告なし

- [ ] **Step 6: コミット**

```bash
git commit -am "Document choosing the relaxation method and add tests with the real OrbMol-v2 model"
```

---

## この計画の外（設計書 6 節）

組み込み後の計測（1 時間あたりの候補数を `n_procs` 4 と 8 で、ピークメモリ、Warp キャッシュの競合、本パイプラインの候補での rejected/failed 率、60 ステップでの順位の保たれ方）は使い捨ての試験として別に行い、結果を `csp_guide.rst` の「計算方法の選び方」「時間の目安」に反映する。
