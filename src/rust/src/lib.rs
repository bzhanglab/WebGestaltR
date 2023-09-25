use extendr_api::prelude::*;
use webgestalt_lib::methods::*;

/// Return string `"Hello world!"` to R.
/// @export
#[extendr]
fn hello_world() -> &'static str {
    "Hello world!"
}

/// Run GSEA using rust library
/// @export
#[extendr]
fn gsea_rust() -> () {
    // webgestalt_lib::methods::gsea::
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod WebGestaltR;
    fn hello_world;
}
