# StarNest examples

What is here, and what is deliberately not.

| File | Purpose |
|---|---|
| `config.example.yaml` | Batch config for `python -m starnest run` (one cluster, UPK 545, as a template) |
| `literature_table.template.csv` | Column layout of the independent literature table `verify.py` needs (header only, no values) |
| `verify.py` | Compares a run's `summary.csv` and member lists with that table, with tolerances fixed at the top of the script |

**No catalogues and no expected outputs are committed.**

- Gaia cone-search CSVs are tens of MB; the repository keeps `data/` out of git. Download each cluster as described in the main README, section 3.0.
- Age, E(B-V), masses, the mass-function slope and the binary-corrected masses depend on the PARSEC isochrone. They are only meaningful from a run with the real grid (`python build_grid.py`, main README section 12.1) or a live CMD request. A run on a mock grid is refused by default and, if forced, is flagged in `summary.csv`.
- Expected outputs are therefore added only from such a run, together with the `run_manifest.json` that produced them (parameters, package versions, input SHA-256), so a reader can see where each number came from.

## Run order for the six study clusters (Haffner 14, Melotte 111, NGC 2910, Ruprecht 147, UBC 52, UPK 545)

1. `pip install -e ".[api]"` and `python build_grid.py` (once, needs internet).
2. Put one cone-search CSV per cluster under `data/` and copy `config.example.yaml` to `my_config.yaml`. Add one entry per cluster under `clusters:`. `cmd.age` (log10 age/yr) is **required** by the batch runner and is an input, not a result: take it from a source you cite, and record that source in your literature table.
3. `python -m starnest run my_config.yaml --out starnest_output`
4. Fill a copy of `literature_table.template.csv`: one row per cluster, from a published catalogue that is **not** StarNest's own output, with the paper or catalogue name in `source` (`verify.py` refuses a row without it). Columns you leave empty are skipped, not passed.
5. `python examples/verify.py starnest_output/summary.csv my_literature_table.csv --csv checks.csv`

`verify.py` exit codes: 0 = nothing failed, 1 = at least one value is outside its tolerance or missing, 2 = unusable input.

What it checks, and why not more: the number of members (N), the distance (d), E(B-V), and the median proper motion and parallax of the members. `logAge` is not checked because it is your own input. The IRSA/SFD E(B-V) is the total dust along the line of sight, so for a nearby cluster it can legitimately exceed the literature value and fail the check; that is a real limitation of the method, not a bug (main README, section 7).
