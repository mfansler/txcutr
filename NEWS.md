# txcutr 0.3.2

NEW FEATURES

* Improved vignette

SIGNIFICANT USER-VISIBLE CHANGES

* None.

BUG FIXES

* None.

# txcutr 0.3.1

NEW FEATURES

* Compressed outputs.
* Tests for proper handling of transitive merging. Overlaps that merge `A -> B` 
  and `B -> C`, but not `A -> C`, will output `A -> C` and `B -> C`. That is, 
  transitivity is applied and the final output will always use the distal most
  transcript in a chain as the final output.

SIGNIFICANT USER-VISIBLE CHANGES

* All `export*()` methods now include automatic detection of `.gz` filenames, 
  which toggles the use of compressed (`gzip`) exports.

BUG FIXES

* None.

# txcutr 0.3.0

NEW FEATURES

* Merge table generation and exporting.

SIGNIFICANT USER-VISIBLE CHANGES

* Adds `generateMergeTable()` and `exportMergeTable()` for creating a merge 
  table for transcripts that are not separated by a thresholded distance. 
  Such files can be used by transcript quantification tools to specify what
  transcripts should be merged.

BUG FIXES

* None.

# txcutr 0.2.2

NEW FEATURES

* None.

SIGNIFICANT USER-VISIBLE CHANGES

* None.

BUG FIXES

* The `BPPARAM` was not being passed through to internal `bplapply` calls.

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
