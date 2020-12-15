from Bio import SeqIO
import functools
import option_helpers as opth
import file_tools


def apply_function_to_sequences(filename, function, seq_type='fasta'):
    with open(filename) as handle:
        for seq_rec in SeqIO.parse(handle, seq_type):
            function(seq_rec)


def apply_function_and_output(input_filename, output_filename, function, seq_type='fasta'):
    with open(output_filename, 'w') as output_handle:
        filter_function = functools.partial(function, output_handle, seq_type)
        apply_function_to_sequences(input_filename, filter_function, seq_type)


class BinarySearcher(object):
    def __init__(self, id_list):
        self.ids = sorted(id_list)

    def binary_search(self, t):
        """
        Take from https://code.activestate.com/recipes/81188-binary-search/
        :param t:
        :return:
        """
        min1 = 0
        max1 = len(self.ids) - 1
        while True:
            if max1 < min1:
                return m
            m = (min1 + max1) // 2
            if self.ids[m] < t:
                min1 = m + 1
            elif self.ids[m] > t:
                max1 = m - 1
            else:
                return m

    def can_find_near_match(self, in_key):
        index = self.binary_search(in_key)
        near_match = False
        if index > 0:
            if in_key.find(self.ids[index - 1]) >= 0:
                near_match = True
        if index < len(self.ids) - 1:
            if in_key.find(self.ids[index + 1]) >= 0:
                near_match = True
        if in_key.find(self.ids[index]) >= 0:
            near_match = True
        return near_match


class CommandFilter(object):
    @staticmethod
    def get_description():
        return 'Filter the sequence to keep (or remove) a list of sequences by id.'

    @staticmethod
    def get_option_data():
        options = opth.default_option_map_input_output()
        options['input']['description'] = 'An input sequence file.'
        options['output']['description'] = 'An output sequence file.'
        options['list'] = {'order': 3,
                           'short': 'l',
                           'long': 'list',
                           'input_name': '<id_list>',
                           'description': "A list of sequence ids in the sequence file.",
                           'optional': False}
        options['remove'] = {'order': 4,
                             'short': 'r',
                             'long': 'remove',
                             'input_name': None,
                             'description': "Remove the sequences in <id_list> rather than keep them.",
                             'optional': True}
        options['near'] = {'order': 5,
                           'short': 'n',
                           'long': 'near',
                           'input_name': None,
                           'description': "Near Match. Select a sequence if it contains the id "
                                          "(appears as a substring).",
                           'optional': True}
        return options

    @staticmethod
    def get_command_data(order=0):
        command_map = {'description': CommandFilter.get_description(),
                       'options': CommandFilter.get_option_data(),
                       'order': order}
        return command_map

    @staticmethod
    def run_program(argument_map, print_usage_function):
        filename = opth.validate_required('input', argument_map, print_usage_function)
        list_filename = opth.validate_required('list', argument_map, print_usage_function)
        output_filename = opth.validate_required('output', argument_map, print_usage_function)
        seq_type = 'fasta'
        if file_tools.looks_fastq(filename):
            seq_type = 'fastq'

        # get the ids from the list file
        id_set = file_tools.parse_lines_to_set(list_filename)

        # If near match is needed, set the flag and create the BinarySearcher
        if not opth.has_option('near', argument_map):
            near = False
        else:
            searcher = BinarySearcher(list(id_set))
            near = True
        # set the remove mode
        remove_mode = opth.has_option('remove', argument_map)

        def is_found(seq, s_type='fasta'):
            if s_type == 'fastq':
                s_id = file_tools.trim_pair_number(seq.id)
            else:
                s_id = seq.id

            if near:
                return searcher.can_find_near_match(s_id)
            else:
                return s_id in id_set

        def decide_and_write_sequence(output_handle, s_type, seq_rec):
            found = is_found(seq_rec, s_type)
            if (found and not remove_mode) or (remove_mode and not found):
                SeqIO.write(seq_rec, output_handle, s_type)
        # pass the function to the application helper
        apply_function_and_output(filename, output_filename, decide_and_write_sequence, seq_type)


class CommandLengthFilter(object):
    @staticmethod
    def get_description():
        return 'Filter the sequence file to keep (or remove) sequence based on their length.'

    @staticmethod
    def get_option_data():
        options = opth.default_option_map_input_output()
        options['input']['description'] = 'An input sequence file.'
        options['output']['description'] = 'An output sequence file.'
        options['length'] = {'order': 3,
                             'short': 'l',
                             'long': 'length',
                             'input_name': '<sequence_length>',
                             'description': "The minimum length to be kept of sequences.",
                             'optional': False}
        options['greater'] = {'order': 4,
                              'short': 'g',
                              'long': 'greater',
                              'input_name': None,
                              'description': "Change filtering to remove sequences greater than <sequence_length>",
                              'optional': True}
        return options

    @staticmethod
    def get_command_data(order=0):
        command_map = {'description': CommandLengthFilter.get_description(),
                       'options': CommandLengthFilter.get_option_data(),
                       'order': order}
        return command_map

    @staticmethod
    def run_program(argument_map, print_usage_function):
        filename = opth.validate_required('input', argument_map, print_usage_function)
        output_filename = opth.validate_required('output', argument_map, print_usage_function)
        length = opth.validate_required_int('length', argument_map, print_usage_function)
        seq_type = 'fasta'
        if file_tools.looks_fastq(filename):
            seq_type = 'fastq'
        filter_greater = opth.has_option('greater', argument_map)

        def keep_sequence(outfile_handle, input_seq_type, seq_rec):
            if (not filter_greater and len(seq_rec) >= length) or (filter_greater and len(seq_rec) < length):
                SeqIO.write(seq_rec, outfile_handle, input_seq_type)
        apply_function_and_output(filename, output_filename, keep_sequence, seq_type)








