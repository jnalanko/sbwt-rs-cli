mod common;
mod build;
mod set_operations;

use clap::Parser;

/// Calls sbwt through the CLI and checks outputs.
#[derive(Parser)]
#[command(name = "sbwt-cli-test-harness", arg_required_else_help = true)]
struct Harness {
    #[command(subcommand)]
    command: Subcommand,
}

#[derive(clap::Subcommand)]
enum Subcommand {
    Build(build::Args),
    SetOperations(set_operations::Args),
}

fn main() {
    match Harness::parse().command {
        Subcommand::Build(args) => build::run(args),
        Subcommand::SetOperations(args) => set_operations::run(args),
    }
}
