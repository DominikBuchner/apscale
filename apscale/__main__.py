import argparse, sys
from apscale import a_create_project, b_pe_merging, c_primer_trimming, d_quality_filtering, e_dereplication_pooling, f_otu_clustering, g_denoising

## main function for the command line interface
def main():
    ## initialize the parse and display default behavior if called without arguments
    formatter = lambda prog: argparse.HelpFormatter(prog, max_help_position = 35)
    parser = argparse.ArgumentParser(prog = 'apscale', description = """Advanced Pipeline for Simple yet Comprehensive AnaLysEs of DNA metabarcoding data,
                                                                     see https://github.com/DominikBuchner/apscale for detailed help.""",
                                                                     formatter_class = formatter)
    # display help when no argument is called
    # set default behavior in case no argument is provided
    parser.set_defaults(func = lambda x: parser.print_help())

    ## add an create project group and the only argument
    create = parser.add_argument_group('Creating a project', description = """Creates a new apscale project in the current working directory""")
    create.add_argument('--create_project', metavar = 'NAME', default = argparse.SUPPRESS, help = 'Creates a new apscale project with the name provided')

    ## add argument groups for the different modules
    modules = parser.add_argument_group('Running a module', description = """Run the apscale pipeline or any specified module. Providing a PATH
                                                                             is optional. If no path is provided apscale will run in the
                                                                             current working directory.""")

    modules.add_argument('--run_apscale', metavar = 'PATH', nargs = '?', const = False, default = argparse.SUPPRESS, help = 'Run the entire pipeline.')
    modules.add_argument('--pe_merging', metavar = 'PATH', nargs = '?', const = False,  default = argparse.SUPPRESS, help = 'Run the pe_merging module.')
    modules.add_argument('--primer_trimming', metavar = 'PATH', nargs = '?', const = False,  default = argparse.SUPPRESS,  help = 'Run the primer_trimimng module.')
    modules.add_argument('--quality_filtering', metavar = 'PATH', nargs = '?', const = False,  default = argparse.SUPPRESS, help = 'Run the quality_filtering module.')
    modules.add_argument('--dereplication_pooling', metavar = 'PATH', nargs = '?', const = False,  default = argparse.SUPPRESS, help = 'Run the dereplication_pooling module.')
    modules.add_argument('--otu_clustering', metavar = 'PATH', nargs = '?', const = False,  default = argparse.SUPPRESS, help = 'Run the otu_clustering module.')
    modules.add_argument('--denoising', metavar = 'PATH', nargs = '?', const = False, default = argparse.SUPPRESS, help = 'Run the denoising module.')

    ## parse command line arguments
    args = parser.parse_args()

    ## main logic of the script
    if 'create_project' in args:
        a_create_project.create_project(args.create_project)

    if 'run_apscale' in args:
        if not args.run_apscale:
            b_pe_merging.main()
            c_primer_trimming.main()
            d_quality_filtering.main()
            e_dereplication_pooling.main()
            f_otu_clustering.main()
            g_denoising.main()
        else:
            b_pe_merging.main(args.run_apscale)
            c_primer_trimming.main(args.run_apscale)
            d_quality_filtering.main(args.run_apscale)
            e_dereplication_pooling.main(args.run_apscale)
            f_otu_clustering.main(args.run_apscale)
            g_denoising.main(args.run_apscale)

    ## check if a module was called, then check if an additional argument was called
    if 'pe_merging' in args:
        if not args.pe_merging:
            b_pe_merging.main()
        else:
            b_pe_merging.main(args.pe_merging)

    ## check if a module was called, then check if an additional argument was called
    if 'primer_trimming' in args:
        if not args.primer_trimming:
            c_primer_trimming.main()
        else:
            c_primer_trimming.main(args.primer_trimming)

    ## check if a module was called, then check if an additional argument was called
    if 'quality_filtering' in args:
        if not args.quality_filtering:
            d_quality_filtering.main()
        else:
            d_quality_filtering.main(args.quality_filtering)

    ## check if a module was called, then check if an additional argument was called
    if 'dereplication_pooling' in args:
        if not args.dereplication_pooling:
            e_dereplication_pooling.main()
        else:
            e_dereplication_pooling.main(args.dereplication_pooling)

    ## check if a module was called, then check if an additional argument was called
    if 'otu_clustering' in args:
        if not args.otu_clustering:
            f_otu_clustering.main()
        else:
            f_otu_clustering.main(args.otu_clustering)

    ## check if a module was called, then check if an additional argument was called
    if 'denoising' in args:
        if not args.denoising:
            g_denoising.main()
        else:
            g_denoising.main(args.denoising)

    ## print help if no argument is provided
    if len(sys.argv) == 1:
        args.func(args)
        sys.exit()

if __name__ == "__main__":
    main()
