import argparse
from apscale import A_create_project

parser = argparse.ArgumentParser(prog = 'apscale', description = 'To be done!')
parser.set_defaults(func = lambda x: parser.print_help())

subparser = parser.add_subparsers(dest = 'test')

parser_create_project = subparser.add_parser('create project')

args = parser.parse_args()
args.func(args)
