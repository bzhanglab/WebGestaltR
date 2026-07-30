/* Keep libc process-termination symbols out of the shared object.
 *
 * Rust's standard library links std::process::{exit,abort} and the runtime cleanup
 * path, whose libc calls appear as undefined `exit`, `_exit` and `abort` in
 * WebGestaltR.so.  No code path in this package reaches them -- ORA, GSEA and NTA
 * never call std::process -- but R CMD check reports the symbols, and their presence
 * alone is what it objects to.
 *
 * The linker's --wrap redirects every reference to __wrap_<name>, so the original
 * names do not appear in the symbol table at all.  The definitions below route any
 * (unreachable) call into R's error mechanism rather than terminating the R session.
 *
 * configure adds the flags only where the linker supports --wrap; see src/Makevars.in.
 */
#include <R.h>
#include <Rinternals.h>

void __wrap_exit(int status);
void __wrap__exit(int status);
void __wrap_abort(void);

void __wrap_exit(int status)  { Rf_error("compiled code attempted exit(%d)", status); }
void __wrap__exit(int status) { Rf_error("compiled code attempted _exit(%d)", status); }
void __wrap_abort(void)       { Rf_error("compiled code attempted abort()"); }
