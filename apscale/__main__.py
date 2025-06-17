import argparse, sys, subprocess, duckdb
from pathlib import Path
from apscale import (
    a_create_project,
    Apscale_analyze,
    b_pe_merging,
    c_primer_trimming,
    d_quality_filtering,
    e_dereplication,
    f_denoising,
    g_swarm_clustering,
    h_replicate_merging,
    i_nc_removal,
    j_generate_read_table,
)


## main function for the command line interface
def main():
    print(duckdb.connect(":memory"))
    ## initialize the parse and display default behavior if called without arguments
    formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position=35)
    parser = argparse.ArgumentParser(
        prog="apscale",
        description="""Advanced Pipeline for Simple yet Comprehensive AnaLysEs of DNA metabarcoding data,
                                                                     see https://github.com/DominikBuchner/apscale for detailed help.""",
        formatter_class=formatter,
    )
    # display help when no argument is called
    # set default behavior in case no argument is provided
    parser.set_defaults(func=lambda x: parser.print_help())

    ## add an create project group and the only argument
    create = parser.add_argument_group(
        "Creating a project",
        description="""Creates a new apscale project in the current working directory""",
    )
    create.add_argument(
        "--create_project",
        metavar="NAME",
        default=argparse.SUPPRESS,
        help="Creates a new apscale project with the name provided",
    )

    ## add argument groups for the different modules
    modules = parser.add_argument_group(
        "Running a module",
        description="""Run the apscale pipeline or any specified module. Providing a PATH
                                                                             is optional. If no path is provided apscale will run in the
                                                                             current working directory.""",
    )

    modules.add_argument(
        "--run_apscale",
        metavar="PATH",
        nargs="?",
        const=False,
        default=argparse.SUPPRESS,
        help="Run the entire pipeline.",
    )
    modules.add_argument(
        "--pe_merging",
        metavar="PATH",
        nargs="?",
        const=False,
        default=argparse.SUPPRESS,
        help="Run the pe_merging module.",
    )
    modules.add_argument(
        "--primer_trimming",
        metavar="PATH",
        nargs="?",
        const=False,
        default=argparse.SUPPRESS,
        help="Run the primer_trimimng module.",
    )
    modules.add_argument(
        "--quality_filtering",
        metavar="PATH",
        nargs="?",
        const=False,
        default=argparse.SUPPRESS,
        help="Run the quality_filtering module.",
    )
    modules.add_argument(
        "--dereplication",
        metavar="PATH",
        nargs="?",
        const=False,
        default=argparse.SUPPRESS,
        help="Run the dereplication_pooling module.",
    )
    modules.add_argument(
        "--denoising",
        metavar="PATH",
        nargs="?",
        const=False,
        default=argparse.SUPPRESS,
        help="Run the denoising module.",
    )

    modules.add_argument(
        "--swarm_clustering",
        metavar="PATH",
        nargs="?",
        const=False,
        default=argparse.SUPPRESS,
        help="Run the swarm clustering module.",
    )

    modules.add_argument(
        "--replicate_merging",
        metavar="PATH",
        nargs="?",
        const=False,
        default=argparse.SUPPRESS,
        help="Run the replicate merging module.",
    )

    modules.add_argument(
        "--nc_removal",
        metavar="PATH",
        nargs="?",
        const=False,
        default=argparse.SUPPRESS,
        help="Run the negative control removal module.",
    )

    modules.add_argument(
        "--generate_read_table",
        metavar="PATH",
        nargs="?",
        const=False,
        default=argparse.SUPPRESS,
        help="Run the read table generation module.",
    )

    modules.add_argument(
        "--analyze",
        metavar="PATH",
        nargs="?",
        const=False,
        default=argparse.SUPPRESS,
        help="Run the analysis module",
    )

    ## parse command line arguments
    args = parser.parse_args()

    ## main logic of the script
    if "create_project" in args:
        a_create_project.create_project(args.create_project)

    if "run_apscale" in args:
        if not args.run_apscale:
            b_pe_merging.main()
            c_primer_trimming.main()
            d_quality_filtering.main()
            e_dereplication.main()
            f_denoising.main()
            g_swarm_clustering.main()
            h_replicate_merging.main()
            i_nc_removal.main()
            j_generate_read_table.main()
        else:
            b_pe_merging.main(args.run_apscale)
            c_primer_trimming.main(args.run_apscale)
            d_quality_filtering.main(args.run_apscale)
            e_dereplication.main(args.run_apscale)
            f_denoising.main(args.run_apscale)
            g_swarm_clustering.main(args.run_apscale)
            h_replicate_merging.main(args.run_apscale)
            i_nc_removal.main(args.run_apscale)
            j_generate_read_table.main(args.run_apscale)

    ## check if a module was called, then check if an additional argument was called
    if "pe_merging" in args:
        if not args.pe_merging:
            b_pe_merging.main()
        else:
            b_pe_merging.main(args.pe_merging)

    ## check if a module was called, then check if an additional argument was called
    if "primer_trimming" in args:
        if not args.primer_trimming:
            c_primer_trimming.main()
        else:
            c_primer_trimming.main(args.primer_trimming)

    ## check if a module was called, then check if an additional argument was called
    if "quality_filtering" in args:
        if not args.quality_filtering:
            d_quality_filtering.main()
        else:
            d_quality_filtering.main(args.quality_filtering)

    ## check if a module was called, then check if an additional argument was called
    if "dereplication" in args:
        if not args.dereplication:
            e_dereplication.main()
        else:
            e_dereplication.main(args.dereplication)

    ## check if a module was called, then check if an additional argument was called
    if "denoising" in args:
        if not args.denoising:
            f_denoising.main()
        else:
            f_denoising.main(args.denoising)

    ## check if a module was called, then check if an additional argument was called
    if "swarm_clustering" in args:
        if not args.swarm_clustering:
            g_swarm_clustering.main()
        else:
            g_swarm_clustering.main(args.swarm_clustering)

    ## check if a module was called, then check if an additional argument was called
    if "replicate_merging" in args:
        if not args.replicate_merging:
            h_replicate_merging.main()
        else:
            h_replicate_merging.main(args.replicate_merging)

    ## check if a module was called, then check if an additional argument was called
    if "nc_removal" in args:
        if not args.nc_removal:
            i_nc_removal.main()
        else:
            i_nc_removal.main(args.nc_removal)

    if "generate_read_table" in args:
        if not args.generate_read_table:
            j_generate_read_table.main()
        else:
            j_generate_read_table.main(args.generate_read_table)

    if "analyze" in args:
        if not args.analyze:
            subprocess.run(
                "streamlit run C:\\Users\\Dominik\\Dokumente\\apscale\\apscale\\Apscale_analyze.py"
            )
        else:
            path = Path(args.analyze)
            subprocess.run(
                f'streamlit run C:\\Users\\Dominik\\Dokumente\\apscale\\apscale\\Apscale_analyze.py -- "{path}"'
            )

    ## print help if no argument is provided
    if len(sys.argv) == 1:
        args.func(args)
        sys.exit()


if __name__ == "__main__":
    main()
