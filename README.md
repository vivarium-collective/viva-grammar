# viva-grammar

<!-- BEGIN:dashboard -->
<!-- END:dashboard -->

A generic, **ABM-agnostic** implementation of the *cell-behavior hypothesis
grammar* — the plain-text modeling language of Johnson, Bergman, …, Fertig &
Macklin, *"Human interpretable grammar encodes multicellular systems biology
models to democratize virtual cell laboratories,"* **Cell** 188:4711–4733 (2025)
— as a [process-bigraph](https://github.com/vivarium-collective/process-bigraph)
module.

Rules read as English (**"In cell type T, signal S increases/decreases behavior
B"**) or as the PhysiCell 8-column CSV; the engine compiles them to the paper's
multivariate Hill-function mathematics (eqs 1–4) and, each step, turns per-cell
**signals** (oxygen, contacts, pressure, …) into per-cell **behaviors** (cycling,
death, migration, type transitions, …). It is designed to plug the *same* rule set
into *different* agent-based models.

> **Clean-room reimplementation, not PhysiCell.** This re-derives the published
> grammar mathematics in portable NumPy; it does not wrap or vendor PhysiCell's
> C++. The reference hypoxia rules and baseline parameters are copied verbatim
> from the paper's `grammar_samples` repo — see
> [`references/PROVENANCE.md`](references/PROVENANCE.md).

## What's here

| Layer | Module | Role |
|---|---|---|
| **Engine** (ABM-agnostic) | `viva_grammar.engine` | signals → behaviors via the Hill math; no ABM imports |
| **IO** | `viva_grammar.io` | PhysiCell rules CSV (canonical) + English-statement parser |
| **Process** | `viva_grammar.steps.GrammarRulesProcess` | the rule engine as a process-bigraph node |
| **Blueprint adapter** | `viva_grammar.steps.BlueprintGrammarProcess` | drop-in for `viva-multicell-blueprint`'s `RulesProcess` |
| **Reference substrate** | `viva_grammar.reference` | a minimal oxygen-field + point-agent ABM for standalone studies |
| **Showcase** | `viva_grammar.composites.hypoxia_motility` | the paper's hypoxia → motility model, end to end |

### Is it a "Step" or a "Process"?

Conceptually it is the "vivarium step that reads rules and executes them." In
process-bigraph terms it is a **`Process`**, not a `Step`, because it needs the
step `interval` to turn physical (per-time) rates into per-step event
probabilities (`Pr = 1 − exp(−rate · interval)`) — and because that makes it a
byte-for-byte drop-in for the blueprint's `RulesProcess` (also a `Process`).

## Install

```bash
# development (editable)
uv venv .venv && source .venv/bin/activate
uv pip install -e ".[dev]"
pytest
```

Once installed, the processes register automatically via
`bigraph_schema.package.discover`; the workspace's `build_core()` also registers
them explicitly for pixi/uv-build editable installs.

## Quick start

```python
from viva_grammar.engine import RuleSet, GrammarEvaluator, SignalContext, CellDefinitions
from viva_grammar.io import parse_rules_csv

rules = RuleSet(parse_rules_csv("tumor,oxygen,decreases,necrosis,0,3.75,8,0"))
ev = GrammarEvaluator(rules, CellDefinitions({"tumor": {"necrosis": {"base": 2.8e-3}}}))

# low oxygen -> necrosis at its base rate; high oxygen -> ~0
print(ev.evaluate(SignalContext("tumor", fields={"oxygen": 0.0}))["necrosis"].value)   # 2.8e-3
print(ev.evaluate(SignalContext("tumor", fields={"oxygen": 100.0}))["necrosis"].value) # ~0
```

Run the showcase study (reference substrate + grammar):

```bash
python studies/hypoxia-motility/sims/run.py
# -> studies/hypoxia-motility/viz/hypoxia_motility.html (interactive)
```

## Showcase: hypoxia → motility

The verbatim PhysiCell hypoxia rules, executed by the generic engine on the
minimal reference substrate, reproduce the paper's **transient post-hypoxic
migration**: a hypoxic core forms, core cells transition to a motile phenotype
and disseminate outward as an invasive plume, the deepest core goes necrotic, and
on reoxygenation the motile cells revert. Acceptance bands (all passing):

| Band | Threshold | Observed |
|---|---|---|
| Phenotype induction | peak motile fraction ≥ 0.15 | ≈ 0.34 |
| Invasive plume | motile radius ≥ 1.2× non-motile (at peak) | ≈ 1.40× |
| Transience / reversion | motile fraction ≤ 0.02 after reoxygenation | → 0.00 |
| Necrotic core | necrotic radius < live radius throughout | central |

> The reference substrate is a deliberately-minimal toy: its oxygen PDE uses a
> reduced diffusivity (5×10³ vs PhysiCell's 10⁵ µm²/min) so a hypoxic core forms
> at the toy's small scale. **The grammar rules and their parameters are
> verbatim**; only the substrate's PDE constants are toy-scaled.

## Plugging into other ABMs — the MISM blueprint

`BlueprintGrammarProcess` matches the `MechanicsBackend` readout/writeback surface
of [`viva-multicell-blueprint`](https://github.com/vivarium-collective/viva-multicell-blueprint)
exactly (inputs `types`/`neighbors`/`shared_surface_by_type`/`field_grid`/…;
outputs `fates`/`scalar`/`divisions`/`deaths`), so the same grammar can drive its
CPM / vertex / center-based backends in place of the built-in `RulesProcess`. It
defaults to `rate_mode="per_step"` to match that engine's `1 − exp(−rate)`
convention.

**Expressiveness boundary:** the grammar combines rules independently per
`(cell_type, behavior)`; it has no rule ordering and no scalar-chaining (unlike
`RulesProcess`'s ordered/`resistance` rules). It is an alternative rule backend
for grammar-expressible models, not a general replacement for ordered rule chains.

## License

MIT.
