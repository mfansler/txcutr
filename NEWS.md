# txcutr 0.2.1

NEW FEATURES

* Added a `NEWS.md` file to track changes to the package.
* Provide more control over parallel execution.

SIGNIFICANT USER-VISIBLE CHANGES

* The `truncateTxome()` method includes an optional `BPPARAM` with which users 
  can pass a specific `BiocParallelParam`. If not provided, it will respect the 
  result of `BiocParallel::bpparam()`, which can be globally set using 
  `BiocParallel::register()`.

BUG FIXES

* None.

# txcutr 0.2.0

NEW FEATURES

* Adds deduplication behavior. Note that deduplication does *not* exclude 
  transcripts from different genes.

SIGNIFICANT USER-VISIBLE CHANGES

* The `truncateTxome()` method now deduplicates transcripts spanning identical 
  ranges after truncation.

BUG FIXES

* None.
