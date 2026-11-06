# Stage 2 — wiring the grammar into viva-multicell-blueprint

Notes for the future worktree session that integrates `viva-grammar` into
[`viva-multicell-blueprint`](https://github.com/vivarium-collective/viva-multicell-blueprint).
The following were **confirmed against blueprint source** by the MISM-blueprint
session (2026-11-06); treat them as settled, not open questions.

## Rate convention — use `per_step` (exact parity)

`RulesProcess.update` is `def update(self, inputs, interval=1)` and fires on
`self._rng.random() < 1.0 - math.exp(-rate)` (`mechanism/rules.py:132`) — it does
**not** multiply `rate` by `interval`. The rate is applied **once per coupled
step**, interval-independent; the runner advances with `interval=1.0`. So
`BlueprintGrammarProcess`'s default `rate_mode="per_step"` (`Pr = 1 − exp(−rate)`)
gives **exact numeric parity**. Do **not** use `per_time` — it coincides only
because `dt == 1` and would silently diverge if a rate were ever interval-scaled.

(The blueprint viz labels the time axis in "days" with `dt_days=1.0`, i.e.
1 step ≈ 1 day — a display convention only; it does not enter the rate math.)

## Composition — strict **swap**, never coexist

There is exactly one rules node: `_shared.py` sets
`RULES_ADDR = "local:RulesProcess"` as the single "rules" slot, wired to write
`fates` / `scalar` / `divisions` / `deaths`. `BlueprintGrammarProcess` occupies
that **same slot** and writes the same four ports. Two writers to those stores
would collide (no merge/arbitration), so the grammar process **replaces**
`RulesProcess`; they never run side by side. `activations` stays untouched — it is
`ImmuneRecruitmentProcess`'s input, a separate node. The death `cause` key is
harmlessly ignored (backends read `cell_id`).

## Boundary — the grammar will NOT reproduce the Sego flu model

`RulesProcess` is "the shared, backend-neutral **ordered** rate-rule chain" and has
a `resistance` rule kind that writes a `scalar` (e.g. `resist`) which later
`rate_expr`s reference. The **Sego influenza blueprint uses** a resistance rule +
field-referencing rate_exprs — i.e. exactly the ordered/resistance chaining the
grammar deliberately does not do (it combines rules independently per
`(cell_type, behavior)`). So:

- The grammar swap is a clean drop-in for **grammar-shaped** blueprints
  (independent per-(type, behavior) rules, no `resist`).
- It will **not** reproduce the current Sego flu rules as-is. **State this
  explicitly in the stage-2 PR** so nobody expects flu-model parity from the swap.
- A good stage-2 demo blueprint is a grammar-shaped one, which also keeps the
  numeric-parity story clean.

## Integration shape

Declare `viva-grammar` under the blueprint's `workspace.yaml` `imports:`
(`package: viva_grammar`, `mode: reference`) — the same pattern used to add
`viva-munk` (center-based paradigm, PR #14). Do all work in a **worktree** of
`viva-multicell-blueprint`; never commit in the shared checkout. Swap
`BlueprintGrammarProcess` into the rules slot for a grammar-shaped study and show
it driving the CPM / vertex / munk backends unchanged.
