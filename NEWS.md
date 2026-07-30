# WebGestaltR 1.0.0

## Breaking changes

- **The default annotation dataset is now WebGestalt 2024.** Previous releases sent no data
  version to the server and therefore received the 2019 dataset. Identical code will now return
  results based on 2024 annotations. Gene sets, category sizes, and consequently p-values and FDRs
  will differ from 0.4.x. This is the intended upgrade, but it is not a silent one: re-run any
  analysis you intend to compare against earlier results.

- **Redundancy reduction is now opt-in, and the default method changed.** 0.4.x always applied
  affinity propagation and weighted set cover. 1.0.0 gates these behind
  `useAffinityPropagation`, `useWeightedSetCover`, and `usekMedoid`, with `usekMedoid = TRUE`
  (`kMedoid_k = 25`) as the default in `WebGestaltR()`. Enrichment statistics are unaffected;
  which categories are reported as representatives will change.

## Bug fixes

- Removed the libc `exit`, `_exit` and `abort` symbols from the compiled shared object. They were
  introduced by Rust's standard library, not by any code path in this package, but `R CMD check`
  reports their presence. Where the linker supports it, references are now redirected so the
  symbols do not appear at all, and any (unreachable) call is raised as an R error rather than
  terminating the R session.

- The summary report now shows the uploaded file name again. A mangled identifier in
  `summaryDescription()` (`interestGeneFilWEBGESTALT_DATA_VERSIONeBase` instead of
  `interestGeneFileBase`) meant the "Interesting list" field rendered empty for all supported
  organisms; only the custom "others" organism was unaffected. Present since 2024-08.

## New features

- Core computations (ORA, GSEA, NTA) are now implemented in Rust for substantially improved
  performance. **This adds a build requirement:** Rust (`cargo`, `rustc` >= 1.63). Binary
  installations from CRAN do not require a local Rust toolchain. On Windows, source installation
  requires the GNU toolchain (`x86_64-pc-windows-gnu`), not MSVC.
- `WebGestaltRMultiOmics()` — multi-omics enrichment across gene and metabolite lists, with
  meta-analysis support.
- `multiswGsea()` — weighted set cover GSEA across multiple lists.
- New arguments to `WebGestaltR()`: `interestGeneNames`, `useWeightedSetCover`,
  `useAffinityPropagation`, `usekMedoid`, `kMedoid_k`, `listName`.

## Compatibility

- No exported functions were removed and no existing argument defaults were changed. Code written
  for 0.4.x will continue to run; see "Breaking changes" above for how its *results* may differ.

# WebGestaltR 0.4.6 (2023-05-31)

- Fixed a bug with R 4.3

# WebGestaltR 0.4.5 (2023-02-10)

- Updated default URLs to HTTPS
- Sort results by enrichment ratio for ties with FDR and P-value
- Fixed a bug in reading GMT with quotes
- Fixed a bug in reading GMT from URL without MIME
- Fixed a bug in GSEA plot path when format is SVG

# WebGestaltR 0.4.4 (2020-07-23)

- Fixed a bug in affinity propagation

# WebGestaltR 0.4.3 (2020-01-16)

Bug fixes. Add a few advanced options for GSEA.

# WebGestaltR 0.4.2 (2019-09-16)

Bug fixes.

# WebGestaltR 0.4.1 (2019-07-03)

Add option to save downloaded data in cache. Bug fixes.

# WebGestaltR 0.4.0 (2019-04-22)

Support of multiple databases for ORA and GSEA. Adjusted column names of the returned data frame.

# WebGestaltR 0.3.1 (2019-03-14)

Bug fixes. The version with NAR update manuscript.

# WebGestaltR 0.3.0 (2019-01-16)

Major updates in HTML report and for 2019 publication

# WebGestaltR 0.1.1

Stable version in WebGestalt 2017 update.

# WebGestaltR 0.0.1

First submission to CRAN.
