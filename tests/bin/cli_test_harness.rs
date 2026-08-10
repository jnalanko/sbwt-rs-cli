mod common;
mod build;
mod set_operations;

fn main() {
    let matches = clap::Command::new("sbwt-cli-test-harness")
        .about("Calls sbwt through the CLI and checks outputs.")
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
