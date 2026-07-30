# cran-comments.md

## Resubmission

This is a second resubmission, addressing the `compiled code` WARNING reported by the Debian
pre-test.

In my previous submission I argued that the finding was a false positive and should be accepted as
such. I understand from your reply that the earlier NOTEs on other Rust-based packages reflected a
gap in the checks rather than acceptable behaviour, and **I withdraw that argument.** The symbols
are now excluded rather than explained — see below.

Also fixed in the earlier resubmission: a scheme-less URI in `README.md`
(`www.webgestalt.org` -> `https://www.webgestalt.org`).

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

* macOS 26.6 (aarch64-apple-darwin23), R 4.6.1 — `R CMD check --as-cran`
* Ubuntu 24.04 (aarch64-unknown-linux-gnu), R 4.5.0, in Docker — `R CMD check --as-cran`,
  used to reproduce the Debian pre-test finding and verify that it is resolved


## Excluding the libc exit/abort symbols

The previous build reported:

```
File 'WebGestaltR/libs/WebGestaltR.so':
  Found '_exit', 'abort', 'exit', possibly from ... (C)
```

**These symbols are now absent from the shared object.**

They originated in Rust's standard library — `std::process::{exit,abort}` and the runtime cleanup
path are linked in as part of the panic runtime, and reference the libc entry points in turn.
Nothing in this package or in its `webgestalt_lib` dependency calls them.

Where the linker supports it, `configure` now adds:

```
-Wl,--wrap=exit -Wl,--wrap=_exit -Wl,--wrap=abort
```

so every reference is redirected to `__wrap_<name>` and the original names do not appear in the
symbol table at all. `src/wrapstubs.c` supplies the three definitions; each raises an R error
rather than terminating the process, so even an unreachable call cannot end the R session.

`configure` probes for `--wrap` support and omits the flags where the linker lacks it, leaving the
build unchanged on such platforms.

Verified with `R CMD check --as-cran`:

* Ubuntu 24.04 (aarch64-unknown-linux-gnu), R 4.5.0 — `checking compiled code ... OK`
* macOS 26.6 (aarch64-apple-darwin23), R 4.6.1 — `checking compiled code ... OK`

I also confirmed with `tools:::check_compiled_code()` directly against the installed package that no
findings remain, and that ORA, GSEA and NTA analyses still run correctly against production data.

## R CMD check results

0 errors | 0 warnings | 1 note.

The note is the maintainer change described above.

Local checks additionally report items that reflect tools absent on the test machines rather than
problems in the package: `checkbashisms` and `V8` on macOS, and `qpdf` on the Linux container.

## System requirements

The package requires Rust (`cargo`, `rustc` >= 1.63) to build from source. Cargo dependencies are
vendored in `src/rust/vendor.tar.xz` and the build runs with `--offline` and a package-local
`CARGO_HOME` unless `NOT_CRAN=true`, so no network access or writes outside the package directory
occur during installation.

On Windows, source installation requires the GNU toolchain (`x86_64-pc-windows-gnu`). This is
documented in the README and in the installation vignette.

## Downstream dependencies

There are no reverse dependencies on CRAN.
