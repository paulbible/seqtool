#!/bin/env python3
from pprint import pprint
import functools
import option_helpers as opth
# command imports
import simple_bio
import seq_filter

def generate_subcommands():
    commands = {'count': simple_bio.CommandCount.get_command_data(order=1),
                'ids': simple_bio.CommandIds.get_command_data(order=2),
                'lengths': simple_bio.CommandLengths.get_command_data(order=3),
                'filter': seq_filter.CommandFilter.get_command_data(order=4),
                'sample': simple_bio.CommandSample.get_command_data(order=5),
                'len-filter': seq_filter.CommandLengthFilter.get_command_data(order=6),
                'edit': simple_bio.CommandEdit.get_command_data(order=7),
                'split': simple_bio.CommandSplit.get_command_data(order=8),
                'prefix': simple_bio.CommandPrefix.get_command_data(order=9),
                'partition': simple_bio.CommandPartition.get_command_data(order=10),
                'base-comp': simple_bio.CommandBaseComposition.get_command_data(order=11),
                'bed-fasta': simple_bio.CommandBedFasta.get_command_data(order=12),
                'to-fasta': simple_bio.CommandToFasta.get_command_data(order=13)}
    return commands


def main():
    subcommand_info = generate_subcommands()
    description = 'A general tool for working with sequences.'
    print_usage_w_subcommands = opth.print_subcommands_usage_maker(subcommand_info, description)

    # extract the command and parse the options
    command = opth.parse_command_options(subcommand_info, print_usage_w_subcommands)
    print_command_usage = functools.partial(print_usage_w_subcommands, command)
    parse_function = opth.parse_options_maker(subcommand_info[command]['options'], print_command_usage, arg_start=2)
    command_options = parse_function()

    if command == 'count':
        simple_bio.CommandCount.run_program(command_options, print_command_usage)
    elif command == 'ids':
        simple_bio.CommandIds.run_program(command_options, print_command_usage)
    elif command == 'lengths':
        simple_bio.CommandLengths.run_program(command_options, print_command_usage)
    elif command == 'filter':
        seq_filter.CommandFilter.run_program(command_options, print_command_usage)
    elif command == 'sample':
        simple_bio.CommandSample.run_program(command_options, print_command_usage)
    elif command == 'len-filter':
        seq_filter.CommandLengthFilter.run_program(command_options, print_command_usage)
    elif command == 'edit':
        simple_bio.CommandEdit.run_program(command_options, print_command_usage)
    elif command == 'split':
        simple_bio.CommandSplit.run_program(command_options, print_command_usage)
    elif command == 'prefix':
        simple_bio.CommandPrefix.run_program(command_options, print_command_usage)
    elif command == 'partition':
        simple_bio.CommandPartition.run_program(command_options, print_command_usage)
    elif command == 'base-comp':
        simple_bio.CommandBaseComposition.run_program(command_options, print_command_usage)
    elif command == 'bed-fasta':
        simple_bio.CommandBedFasta.run_program(command_options, print_command_usage)
    elif command == 'to-fasta':
        simple_bio.CommandToFasta.run_program(command_options, print_command_usage)


if __name__ == '__main__':
    main()
