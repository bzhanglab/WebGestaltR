# cran-comments.md

## Submission summary

This is a major release (0.4.6 -> 1.0.0). Core computations (ORA, GSEA, NTA) are now implemented
in Rust for improved performance, and the package gains multi-omics support
(`WebGestaltRMultiOmics()`).

Two user-visible behaviour changes are documented prominently in `NEWS.md`:

* The default annotation dataset is now WebGestalt 2024. Earlier releases sent no data version to
  the server and therefore received the 2019 dataset, so identical code will return results based
  on newer annotations.
* Redundancy reduction is now opt-in, with k-medoid as the default method rather than affinity
  propagation plus weighted set cover.

No exported functions were removed and no existing argument defaults were changed, so code written
for 0.4.x continues to run.

## Maintainer change

This release changes the maintainer, which produces the expected NOTE.

The previous maintainer, Yuxing Liao, and the author of the 1.0.0 development work, John
Elizarraras, have both left Baylor College of Medicine. I am a co-author of the package — listed
in `Authors@R` in the current CRAN release (0.4.6) — and I now maintain WebGestalt at the Zhang
lab. The package is developed under the lab's GitHub organisation at
<https://github.com/bzhanglab/WebGestaltR>.

I have contacted the previous maintainer regarding this transfer. The lab's principal
investigator, Bing Zhang (a co-author on the WebGestalt 2024 publication cited in `inst/CITATION`),
can also confirm it if that is helpful.

The maintainer address is a personal one rather than an institutional one, deliberately, so that
CRAN correspondence continues to reach me.

## Test environments

* local: macOS 26.6 (aarch64-apple-darwin23), R 4.6.1 — `R CMD check --as-cran`

## R CMD check results

0 errors | 0 warnings | 1 note

The note is the maintainer change described above.

Two further items appear in the local check but reflect tools that are absent on the test machine
rather than problems in the package: `checkbashisms` (not available on macOS) and `V8` (used only
for math-rendering checks of the HTML manual). Both should run normally on the CRAN builders.

## System requirements

The package requires Rust (`cargo`, `rustc` >= 1.63) to build from source. Cargo dependencies are
vendored in `src/rust/vendor.tar.xz` and the build runs with `--offline` and a package-local
`CARGO_HOME` unless `NOT_CRAN=true`, so no network access or writes outside the package directory
occur during installation.

On Windows, source installation requires the GNU toolchain (`x86_64-pc-windows-gnu`). This is
documented in the README and in the installation vignette.

## Downstream dependencies

There are no reverse dependencies on CRAN.
