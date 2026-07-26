# Trace schedulers in HEBCPF v5

The asynchronous `parfeval` drivers in both v5 implementations expose three
canonical trace schedulers: `scan`, `bandit`, and `novelty`. `bandit` is the
default. The legacy internal name `diverse` remains a deprecated alias for
`bandit` so old scripts continue to run.

```matlab
global HEBCPOLICY
HEBCPOLICY = 'bandit';        % default; alternatively scan or novelty
[result,solutions] = run_merged_case('case57','parfeval');
```

`trace_policy_config` normalizes the name in `main.m`, wrappers, batch scripts,
and queue drivers. Unknown names fail immediately instead of silently selecting
a different policy.

## Shared queue and invariants

The search state is the same under every scheduler. `VBook(s,e)` records
whether the continuation curve containing solution `s` has been covered on
equation `e`. A pair is eligible when that entry is zero and equation `e` has
no in-flight task. At most one trace per equation is in flight.

Each equation owns a cursor into solution indices. Selection advances the
cursor past already-covered solutions; newly discovered solutions are reached
later without appending them to per-equation pending lists. This removes the
v4 auxiliary `O(ns*neq)` pending-list duplication. `VBook` itself necessarily
remains `O(ns*neq)` because it is the coverage matrix.

All policies use the same continuation kernels, solution deduplication,
coverage stamping, termination test, and checkpoint format. Thus the intended
difference is dispatch order and the resulting anytime profile. Parallel
completion order is nondeterministic, so trace-to-target and wall time are
measurements rather than invariants.

## `scan`: v4-compatible rotating fairness

`scan` visits equations from a rotating pointer and selects the first non-busy
equation with an uncovered solution. It is the direct v4 scheduling baseline
and is the appropriate choice for reproducing v4 queue order as closely as
parallel completion permits.

Cursor checks are constant work individually and amortize across skipped
covered pairs; auxiliary scheduler state is `O(neq)`. `scan` adds no gain sort
or geometric novelty calculation.

## `bandit`: learned equation productivity (default)

`bandit` learns which equations have recently exposed new solutions. Each
equation has a discounted gain

```text
g[e] <- 0.7*g[e] + 0.3*(new solutions from the completed trace on e)
```

Gains start optimistically at 5 so equations are explored. For each dispatch,
equations are sorted by descending gain and the first eligible equation is
used; its solution is obtained through the shared cursor. Sorting costs
`O(neq log neq)` per dispatch and the gain vector costs `O(neq)` memory.

After `2*neq` consecutive completed traces discover nothing, the learned order
has no useful signal for the mechanical coverage tail. The driver therefore
hands off one-way to rotating `scan` order for the remainder of that run.
Discovery resets the stall counter.

## `novelty`: productivity plus geometric separation

`novelty` uses the same gain-ranked equation choice as `bandit`, but changes
the starting solution within the selected equation. It samples up to 96
eligible candidates and up to 96 already-covered references, projects them to
eight dimensions with a fixed random projection, and chooses the candidate
farthest from the covered set (including antipodal references).

The sampling work is capped at 400 draws per dispatch. With the configured
caps, the distance calculation is bounded by
`O(k_candidate*k_reference*k_projection)` rather than scaling over the whole
solution set. The projection costs `O(numvar*k_projection)` memory; candidate
and reference buffers are bounded. Dedicated fixed-seed streams keep novelty
sampling independent of unrelated MATLAB random-number use, although worker
completion order still makes a parallel run nondeterministic.

`novelty` is intended to test whether geometrically separated starts improve
early discovery. Its extra projection, sampling, and distance work can outweigh
that benefit, particularly on small cases or during the coverage tail.

## Choosing a scheduler

- Use `bandit` for the v5 default balance of early discovery and bounded
  scheduler overhead.
- Use `scan` for direct comparison with the v4 queue policy and for the lowest
  selection overhead.
- Use `novelty` when early geometric diversity matters enough to justify the
  additional bounded selection work.

The full recommendation must be based on the measured case and hardware. The
release benchmark reports average wall time across three repeats, traces and
wall time to 90% of the final solution count, total traces, solution counts,
worker trace-time statistics, scheduler selection time, and solution-set
agreement for all bundled cases through case57.

## Measured v5 performance

The release benchmark completed 180/180 timing runs on MATLAB R2022a, Windows
x64, with 23 workers. Summed per-case mean exhaustive wall time was nearly
identical: 606.947 s for `scan`, 607.184 s for `bandit`, and 606.698 s for
`novelty`. The scheduler distinction is therefore mainly an anytime result in
this suite:

| Policy | Sum traces to 90% | Change from scan | Selection time | Selection share |
| --- | ---: | ---: | ---: | ---: |
| `scan` | 35,975 | baseline | 3.573 s | 0.59% |
| `bandit` | 16,611 | 53.8% lower | 0.980 s | 0.16% |
| `novelty` | 11,975 | 66.7% lower | 34.163 s | 5.63% |

On case57, mean traces to 90% were 19,481 (`scan`), 6,996 (`bandit`), and
4,594 (`novelty`); exhaustive wall times remained 333--337 s. On case57mod,
the corresponding trace counts were 9,203, 4,425, and 2,829. This supports
`bandit` as the default: it captures most of the early-discovery benefit with
far less selection work than `novelty`. `novelty` is the stronger experimental
choice when the earliest broad discovery matters more than selection cost.

All policies returned the expected solution counts. A separate 40-run
order-independent validation found a worst symmetric nearest-solution distance
of `2.66e-9`, below the `4e-7` tolerance. Total trace count varied slightly in
some cases because asynchronous work can finish after another task makes it
redundant; it is consequently reported as a measured outcome.

See `SCHEDULER_BENCHMARK_v5.md` for all cases and plots.

## Checkpoint behavior

Checkpoints store the common solver state, trace history, normalized policy
name, and `eq_gain_d`. Any canonical scheduler can resume another scheduler's
checkpoint. `bandit` and `novelty` warm-start compatible learned gains; `scan`
ignores them. Policy-local equation cursors, busy flags, stall counters, random
streams, and novelty projections are rebuilt from the saved coverage state.

`test_checkpoint_policy_compatibility_v5.m` checks scan→bandit,
bandit→novelty, and novelty→scan transitions under both MEX and pure MATLAB.
The deterministic test controls are:

- `HEBCPF_CHECKPOINT_TRACE_INTERVAL=<positive integer>`
- `HEBCPF_STOP_AFTER_CHECKPOINT=true`

Leave both unset for normal runs; the original time-based production cadence
then remains active.

The scheduler applies only to `runVBook_hybrid_parfeval.m`. The serial,
`parfor`, and retained row-barrier drivers use their native traversal orders.
