import argparse
from apscale import A_create_project

## initialize the parse and display default behavior if called without arguments
parser = argparse.ArgumentParser(prog = 'apscale', description = 'To be done!')
parser.set_defaults(func = lambda x: parser.print_help())

## parse command line arguments
args = parser.parse_args()

## display help when no argument is called
args.func(args)
