//! Standalone harness (not run by `cargo test`) for exercising the `sbwt` CLI end to end,
//! rather than calling the library builders/set operations directly.
//!
//! Subcommands:
//!   build            drives all four construction algorithms on the same input and checks
//!                    that they produce byte-identical output.
//!   set-operations   builds one SBWT per input file and checks merge/intersect/difference
//!                    between them.
//!
//! Run `cargo run --release --features libsais --bin sbwt-construction-test-harness -- <subcommand> --help`
//! for each subcommand's options.

mod common;
mod build;
mod set_operations;

fn main() {
    let matches = clap::Command::new("sbwt-construction-test-harness")
        .about("Drives the sbwt CLI through construction and set-operation workflows and \
                checks its output.")
        .subcommand_required(true)
        .arg_required_else_help(true)
        .subcommand(build::subcommand())
        .subcommand(set_operations::subcommand())
        .get_matches();

    match matches.subcommand() {
        Some(("build", sub_matches)) => build::run(sub_matches),
        Some(("set-operations", sub_matches)) => set_operations::run(sub_matches),
        _ => unreachable!("clap subcommand_required guarantees one of the above matched"),
    }
}
