# Analysis artifacts

No large or full-study outputs are retained in this directory. No full category matrix, raw simulation CSV collection, plot set, or default `10000 × 500` study result is present.

The retained bounded validation inputs and results are:

- [`../data/reference/publication_cases.csv`](../data/reference/publication_cases.csv) — approved Patient room and Museum bounded regression rows
- [`../docs/VALIDATION.md`](../docs/VALIDATION.md) — status, exact bounded results, counters, provenance, and limitations

The evidence is bounded validation only. Patient room uses `N=10000`, `R=10`; Museum uses `N=10000`, `R=5`. The full default workload has not been run and remains pending explicit authorization.

No analysis scripts or raw CSV files are required to reproduce the retained checks. The end-to-end test invokes the production executable, generates temporary summary and replicate output where required, verifies the results, and removes the temporary output tree.

After configuring and building the `ecai_cpp2_end_to_end_tests` target, run the focused check from Windows with:

```bash
cd "$CPP2"
cmake --build build --target ecai_cpp2_end_to_end_tests
./build/ecai_cpp2_end_to_end_tests
```

A passing focused run reports 52 deterministic regression, wiring, documentation, bounded-comparator, accounting, and cleanup checks. These checks confirm the code behaves as implemented against the retained fixtures described in [`docs/VALIDATION.md`](../docs/VALIDATION.md).
