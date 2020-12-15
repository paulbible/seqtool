# Funtion helpers for working with options
import getopt
import sys
import os
import pprint


def help_option_map():
    '''
    An example of an option. All options should be a map having these fields
    :return:
    '''
    return {'help': {'order': 0,
                     'short': 'h',
                     'long' : 'help',
                     'input_name': None,
                     'description': 'Print this help message',
                     'optional': False}}


def sys_exit(code=0):
    sys.exit(code)


def print_example_option():
    pprint.pprint(help_option_map())


def print_example_subcommand_option():
    options = {'command1': {'description': 'A description of command 1.',
                            'options': opth.default_option_map_input(),
                            'order': 1},
               'command2': {'description': 'A description of command 2.',
                            'options': opth.help_option_map(),
                            'order': 2}}
    pprint.pprint(options)


def default_option_map_input():
    option_map = help_option_map()
    option_map['input'] = {'order': 1,
                           'short': 'i',
                           'long': 'input',
                           'input_name': '<input_file>',
                           'description': 'An input file. [Change Me]',
                           'optional': False}
    return option_map


def default_option_map_input_output():
    option_map = default_option_map_input()
    option_map['output'] = {'order': 2,
                            'short': 'o',
                            'long': 'output',
                            'input_name': '<ouput_file>',
                            'description': 'An output file. [Change Me]',
                            'optional': False}
    return option_map


def opt_to_options_transform(option_map):
    new_map = {}
    for k in option_map:
        opt_rec = option_map[k]
        short_key = '-' + opt_rec['short']
        new_map[short_key] = opt_rec
        long_key = '--' + opt_rec['long']
        new_map[long_key] = opt_rec
    return new_map


def print_usage_maker_base(command_description, options, command_name, epilogue=None):
    command_string = 'usage: %s' % command_name
    command_string += ' '
    option_recs = list(options.values())
    option_recs.sort(key=lambda x: x['order'])
    for opt in option_recs:
        if opt['short'] == 'h':
            continue
        if opt['optional']:
            open_optional = ' ['
            close_optional = ']'
        else:
            open_optional = ''
            close_optional = ''
        if opt['input_name']:
            command_string += open_optional + ' -' + opt['short'] + ' ' + opt['input_name'] + close_optional
        else:
            command_string += open_optional + ' -' + opt['short'] + close_optional
    option_descriptions = 'Options:\n'
    spacer = '    '
    for opt in option_recs:
        opt_str = spacer
        if opt['input_name']:
            opt_str += '-' + opt['short'] + ', --' + opt['long'] + '= ' + opt['input_name']
        else:
            opt_str += '-' + opt['short'] + ', --' + opt['long']
        total_len = len(opt_str)
        if total_len < 50:
            diff = 50 - total_len
            between_spacer = ""
            for i in range(diff):
                between_spacer += ' '
            opt_str += between_spacer
        else:
            # add 50 spaces
            opt_str += '\n                                                  '
        if opt['optional']:
            opt_str += '(Optional) ' + opt['description']
        else:
            opt_str += opt['description']
        option_descriptions += opt_str + '\n'

    def return_func(msg=''):
        if len(msg) != 0:
            print(msg)
        print(command_string)
        print(spacer + command_description)
        print(option_descriptions)
        if epilogue:
            print(epilogue)
    return return_func


def print_usage_maker(command_description, options, epilogue=None):
    return print_usage_maker_base(command_description, options, os.path.basename(sys.argv[0]), epilogue=epilogue)


def print_subcommands_usage_maker(subcommand_options, description):
    command_list = []
    for command in subcommand_options:
        command_list.append((subcommand_options[command]['order'], command))
    command_list.sort()
    command_order = list(map(lambda x: x[1], command_list))

    program_name = os.path.basename(sys.argv[0])
    command_string = 'usage: %s command [options]' % program_name
    command_descriptions = 'Commands:\n'
    spacer = '    '
    for command in command_order:
        opt_str = spacer + command + ':'
        total_length = len(opt_str)
        if total_length < 25:
            diff = 25 - total_length
            between_spacer = ""
            for i in range(diff):
                between_spacer += ' '
            opt_str += between_spacer
        opt_str += subcommand_options[command]['description']
        command_descriptions += opt_str + '\n'

    print_map = {}
    for command in command_order:
        if 'epilogue' in subcommand_options[command]:
            epilogue = subcommand_options[command]['epilogue']
        else:
            epilogue = None
        print_map[command] = print_usage_maker_base(subcommand_options[command]['description'],
                                                    subcommand_options[command]['options'],
                                                    program_name + ' ' + command,
                                                    epilogue=epilogue)

    def print_usage_func(command=None, msg=''):
        if len(msg) != 0:
            print(msg)

        if not command or command not in subcommand_options:
            print(command_string)
            print(spacer + description)
            print(command_descriptions)
            print('run \'%s command -h\' for more information on a given command.' % program_name)
        elif command in subcommand_options:
            print_map[command]()

    return print_usage_func


def parse_options_maker(options, print_usage_func, arg_start=1):
    short_opt_codes = ''
    long_opt_codes = []
    option_recs = list(options.values())
    option_recs.sort(key=lambda x: x['order'])
    for opt_rec in option_recs:
        if opt_rec['input_name']:
            short_opt_codes += opt_rec['short'] + ':'
            long_opt_codes.append(opt_rec['long'] + '=')
        else:
            short_opt_codes += opt_rec['short']
            long_opt_codes.append(opt_rec['long'])
    flag_option_map = opt_to_options_transform(options)

    def parse_func():
        try:
            opts, args = getopt.getopt(sys.argv[arg_start:], short_opt_codes, long_opt_codes)
        except getopt.GetoptError as err:
            print_usage_func(str(err))
            sys.exit(2)
        output_map = {}
        if len(opts) == 0:
            print_usage_func()
            sys.exit(1)
        for opt, arg in opts:
            if opt in ("-h", "-help", "--help"):
                print_usage_func()
                sys.exit(1)
            elif opt in flag_option_map:
                option_name = flag_option_map[opt]['long']
                output_map[option_name] = arg
            else:
                print_usage_func("Error: Option '%s' is not recognized." % opt)
                sys.exit(1)
        for option in options:
            if option == 'help':
                continue
            if not options[option]['optional'] and option not in output_map:
                print_usage_func("Error: Required argument '%s' is missing." % option)
                sys.exit(1)
        return output_map
    return parse_func


def parse_command_options(command_options, print_command_usage_funcs):
    if len(sys.argv) < 2:
        print_command_usage_funcs()
        sys.exit(1)

    command = sys.argv[1]
    if command not in command_options:
        print_command_usage_funcs(msg='ERROR: \'%s\' is not a recognized command' % command)
        sys.exit(1)
    else:
        return command


def validate_required(argument_key, argument_map, print_usage_function):
    if argument_key not in argument_map or not argument_map[argument_key]:
        print_usage_function("Error: The %s argument is required." % argument_key)
        sys.exit(1)
    else:
        return argument_map[argument_key]


def has_option(option, option_map):
    return option in option_map


def is_valid_int(argument_key, argument_map):
    try:
        int(argument_map[argument_key])
    except ValueError:
        return False
    return True


def is_valid_float(argument_key, argument_map):
    try:
        float(argument_map[argument_key])
    except ValueError:
        return False
    return True


def with_default(argument_key, argument_map, default_value):
    if has_option(argument_key, argument_map):
        return argument_map[argument_key]
    else:
        return default_value


def with_default_int(argument_key, argument_map, default_value):
    if has_option(argument_key, argument_map) and is_valid_int(argument_key, argument_map):
        return int(argument_map[argument_key])
    else:
        return default_value


def with_default_float(argument_key, argument_map, default_value):
    if has_option(argument_key, argument_map) and is_valid_float(argument_key, argument_map):
        return float(argument_map[argument_key])
    else:
        return default_value


def validate_required_float(argument_key, argument_map, print_usage_function):
    value = validate_required(argument_key, argument_map, print_usage_function)
    try:
        float(value)
    except ValueError:
        print_usage_function('ERROR: %s does not appear to be a float.' % value)
        sys.exit(1)
    return float(value)


def validate_required_int(argument_key, argument_map, print_usage_function):
    value = validate_required(argument_key, argument_map, print_usage_function)
    try:
        int(value)
    except ValueError:
        print_usage_function('ERROR: %s does not appear to be an integer.' % value)
        sys.exit(1)
    return int(value)
