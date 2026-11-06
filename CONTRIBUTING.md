# Contributing to viva-grammar

## Development setup

uv is required. Install with `brew install uv` or `pip install uv`.

    uv venv .venv
    source .venv/bin/activate
    uv pip install -e ".[dev]"
    pytest

## Scope

`viva-grammar` is a clean-room reimplementation of the cell-behavior hypothesis
grammar (Johnson et al., *Cell* 2025). The engine (`viva_grammar/engine/`) must
stay **ABM-agnostic** — it imports no specific agent-based model. All coupling to
a host ABM lives in `viva_grammar/steps/adapters/`.

## Releasing to PyPI

Tag a commit with `git tag v<VERSION>` and push the tag; `release.yml` publishes
via PyPI trusted publishing. See https://docs.pypi.org/trusted-publishers/.
