## Submission

WebGestaltR 1.0.1 fixes the installation failure reported on the macOS check flavors
and under Additional issues (M1mac).

This arrives one day after 1.0.0 because it corrects check failures rather than adding
features.

## The failure

All four macOS flavors, and the M1mac additional check, failed to install:

```
make: rustc: No such file or directory
make: *** [rust/target/release/libWebGestaltR.a] Error 1
ERROR: compilation failed for package 'WebGestaltR'
```

`src/Makevars.in` ran `rustc --version` as a diagnostic before the recipe extended
`PATH` with `~/.cargo/bin`. Because each line of a make recipe runs in its own shell,
that extension — written on the `cargo build` line — never applied to the version
check. On any system where the Rust toolchain is not on the default `PATH`, the build
aborted before cargo was invoked.

Both `src/Makevars.in` and `src/Makevars.win.in` now define the extended `PATH` once and
use it for every step, including the version check. `Makevars.win.in` had the same
ordering and is corrected as well; it had not failed because the Windows builder happens
to have `rustc` on `PATH`.

## Testing

Reproduced and verified in both directions, with `rustc` present in `~/.cargo/bin` but
absent from `PATH`:

* Debian, R 4.3.2 — fails with the error above before the change; installs after it.
* macOS 26 (arm64), R 4.6 — full `R CMD INSTALL` and `R CMD check --as-cran` pass after
  the change. `checking compiled code` is OK.

The macOS installations had never previously completed, so this is the first
confirmation that omitting the `--wrap` linker flags there is correct: Apple's `ld` has
no `--wrap`, `configure` detects that and emits no flags, and the resulting shared
object passes `checking compiled code` on macOS.

## Additional issues (M1mac)

Same cause, confirmed from
<https://www.stats.ox.ac.uk/pub/bdr/M1mac/WebGestaltR.log>:

```
make[1]: rustc: No such file or directory
```

The log shows cargo being found and its version accepted, and `entrypoint.c` and
`wrapstubs.c` compiling, before the build stops at `rust/target/release/libWebGestaltR.a`.
That is precisely the failure described above: `configure` locates cargo, and then the
`rustc --version` line in `src/Makevars` runs in a shell without `~/.cargo/bin` on
`PATH`. This release fixes it.

## R CMD check results

`R CMD check --as-cran` on macOS 26 (arm64), R 4.6: 0 errors | 1 warning | 1 note.

* WARNING: `A complete check needs the 'checkbashisms' script` — not installed on my
  machine. `configure` and `configure.win` pass `sh -n`.
* NOTE: `Days since last update: 1` — this release exists to fix the reported check
  failures.

Debian, R 4.3.2: installs and loads cleanly.

One NOTE remains on r-devel-windows-x86_64, unchanged by this release:

```
Error in ccE(lines, flags = new_flags, include = include) : 'cc' is not on the path
Calls: <Anonymous> ... lapply -> FUN -> lapply -> FUN -> getFunsHdr -> ccE
```

This is the `checking compiled code` step failing to run rather than reporting a
finding: it looks for a compiler named `cc`, and that builder provides `gcc.exe`. Every
other check on that flavor is OK. I have not been able to reproduce it without an
R-devel Windows machine — please let me know if you would like it addressed and can
suggest what the package should change.
