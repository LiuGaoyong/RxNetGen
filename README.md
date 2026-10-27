# RxNetGen

The reaction network generator & event table exploration.

RxNetGen builds reaction networks and explores event tables for KMC
simulation.  The chemical core (systems, graphs, geometry, dataclasses and
array-api helpers) is provided by the
[`graphatoms`](https://github.com/LiuGaoyong/GraphAtoms) package; RxNetGen
depends on it and provides the application layer: configuration, reaction
network storage/commands, runners, parallel executors, reaction/event classes
and the CLI entry points.

## Installation

```bash
pip install -e .
```

Or with `uv` / `pixi`:

```bash
uv sync
pixi install
```

## CLI entry points

| Command | Description |
| --- | --- |
| `rxnetgen-config` | Resolve and print the Hydra/OmegaConf run configuration |
| `rxnetgen-run` | Launch a run (`run_type=otfkmc` or `run_type=rxngen`) |
| `rxnetgen-network` | Inspect and visualize a stored reaction network |

## Module structure

```
src/rxnetgen/
├── config/     # Hydra/OmegaConf configuration (atoms, bonds, calculator)
├── network/    # Reaction network: scheduler, recorder, metadata, commands
├── parallel/   # Pluggable executors (serial, multiprocessing, joblib, ray, dask, executorlib)
├── reaction/   # Reaction classes and KMC events
├── runner/     # Runners (otfkmc, rxngen) and helpers
└── utils/      # Vendored copy of graphatoms.utils (app-level helpers)
```

> **Note**: `rxnetgen.utils` is a vendored copy of `graphatoms.utils` kept so
> that the application layer can evolve independently.  Changes made upstream
> in `graphatoms.utils` must be synchronized here when they affect the
> migrated modules.

## Running tests

```bash
pytest src/tests -v
```

The tests use the datasets under `src/tests-datasets/`.
