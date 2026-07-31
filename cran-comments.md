## Submission

This release fixes the installation failure reported on the four macOS flavors and
under Additional issues (M1mac):

```
make: rustc: No such file or directory
```

`src/Makevars.in` ran `rustc --version` before the recipe added `~/.cargo/bin` to
`PATH`. Each line of a make recipe runs in its own shell, so that addition — written on
the `cargo build` line — never applied to the version check, and the build stopped
before cargo ran. Both Makevars files now set `PATH` once and use it for every step.

It arrives a day after 1.0.0 because it corrects the reported check failures.

## Test environments

* macOS 26 (arm64), R 4.6 — `R CMD check --as-cran`
* Debian, R 4.3.2 — install and load

Both were run with `rustc` present in `~/.cargo/bin` but absent from `PATH`, which
reproduces the reported failure before the change and passes after it.

## R CMD check results

0 errors | 1 warning | 1 note

* WARNING: `checkbashisms` is not installed here; `configure` and `configure.win` pass
  `sh -n`.
* NOTE: days since last update — see above.

## Note on r-devel-windows-x86_64

The existing NOTE there is `checking compiled code` failing to run rather than
reporting a finding:

```
Error in ccE(...) : 'cc' is not on the path
```

That builder provides `gcc.exe`. Every other check on the flavor is OK, and I could not
reproduce this without an R-devel Windows machine. Happy to make changes if you can
suggest what the package should do differently.
