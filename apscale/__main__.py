import argparse, sys
from apscale import a_create_project

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
    if args.create_project:
        a_create_project.create_project(args.create_project)

    # ## print help if no argument is provided
    # if len(sys.argv) == 1:
    #     args.func(args)
    #     sys.exit()

main()
