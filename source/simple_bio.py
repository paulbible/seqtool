from Bio import Seq, SeqIO, SeqRecord
import option_helpers as opth
import functools
import file_tools
from collections import defaultdict
from random import random, randrange


def apply_function_to_sequences(filename, function, seq_type='fasta'):
    with open(filename) as handle:
        for seq_rec in SeqIO.parse(handle, seq_type):
            function(seq_rec)


def apply_function_and_output(input_filename, output_filename, function, seq_type='fasta'):
    with open(output_filename, 'w') as output_handle:
        filter_function = functools.partial(function, output_handle, seq_type)
        apply_function_to_sequences(input_filename, filter_function, seq_type)


# !---------- template command ----------!
class CommandTemplate(object):
    @staticmethod
    def get_description():
        pass

    @staticmethod
    def get_option_data():
        pass

    @staticmethod
    def get_command_data(order=0):
        pass

    @staticmethod
    def run_program(argument_map, print_usage_function):
        pass


# !---------- 'Count' command ----------!
class CommandCount(object):
    @staticmethod
    def get_description():
        return 'Count the number of sequences in a FASTA file.'

    @staticmethod
    def get_option_data():
        options = opth.default_option_map_input()
        options['input']['description'] = 'An input fasta or fastq file'
        options['input']['input_name'] = '<sequence_file>'
        return options

    @staticmethod
    def get_command_data(order=0):
        command_map = {'description': CommandCount.get_description(),
                       'options': CommandCount.get_option_data(),
                       'order': order}
        return command_map

    @staticmethod
    def run_program(argument_map, print_usage_function):
        opth.validate_required('input', argument_map, print_usage_function)
        filename = argument_map['input']

        seq_type = 'fasta'
        if file_tools.looks_fastq(filename):
            seq_type = 'fastq'

        count = 0
        with open(filename) as handle:
            for seq_rec in SeqIO.parse(handle, seq_type):
                count += 1
        print(count)


# !---------- 'ids' command ----------!
class CommandIds(object):
    @staticmethod
    def get_description():
        return 'Output all sequence ids in a fasta or fastq file.'

    @staticmethod
    def get_option_data():
        options = opth.default_option_map_input()
        options['input']['description'] = 'An input fasta file'
        options['input']['input_name'] = '<sequence_file>'
        return options

    @staticmethod
    def get_command_data(order=0):
        command_map = {'description': CommandIds.get_description(),
                       'options': CommandIds.get_option_data(),
                       'order': order}
        return command_map

    @staticmethod
    def run_program(argument_map, print_usage_function):
        opth.validate_required('input', argument_map, print_usage_function)
        filename = argument_map['input']

        seq_type = 'fasta'
        if file_tools.looks_fastq(filename):
            seq_type = 'fastq'

        def print_func(seq_rec):
            print(seq_rec.id)
        apply_function_to_sequences(filename, print_func, seq_type)


# !---------- 'lengths' command ----------!
class CommandLengths(object):
    @staticmethod
    def get_description():
        return 'Output the sequence id and length for each sequence.'

    @staticmethod
    def get_option_data():
        options = opth.default_option_map_input()
        options['input']['description'] = 'An input fasta file'
        options['input']['input_name'] = '<sequence_file>'
        options['delimiter'] = {'order': 2,
                                'short': 'd',
                                'long': 'delimiter',
                                'input_name': '<delimiter_character>',
                                'description': "A character to use as a delimiter e.g. '-d ,'",
                                'optional': True}
        return options

    @staticmethod
    def get_command_data(order=0):
        command_map = {'description': CommandLengths.get_description(),
                       'options': CommandLengths.get_option_data(),
                       'order': order}
        return command_map

    @staticmethod
    def run_program(argument_map, print_usage_function):
        opth.validate_required('input', argument_map, print_usage_function)
        filename = argument_map['input']

        seq_type = 'fasta'
        if file_tools.looks_fastq(filename):
            seq_type = 'fastq'

        if 'delimiter' not in argument_map or not argument_map['delimiter']:
            delimiter = '\t'
        else:
            delimiter = argument_map['delimiter']

        def print_func(seq_rec):
            print(delimiter.join([seq_rec.id, str(len(seq_rec))]))
        apply_function_to_sequences(filename, print_func, seq_type)


# !---------- 'sample' command ----------!
class CommandSample(object):
    @staticmethod
    def get_description():
        return 'Randomly extract a percentage of sequences from a file.'

    @staticmethod
    def get_option_data():
        options = opth.default_option_map_input_output()
        options['input']['description'] = 'An input sequence file'
        options['output']['description'] = 'An output sequence file'
        options['percentage'] = {'order': 3,
                                 'short': 'p',
                                 'long': 'percentage',
                                 'input_name': '<float_percent>',
                                 'description': 'Add a sequence to the output file if a random value <= p',
                                 'optional': False}
        return options

    @staticmethod
    def get_command_data(order=0):
        command_map = {'description': CommandSample.get_description(),
                       'options': CommandSample.get_option_data(),
                       'order': order}
        return command_map

    @staticmethod
    def run_program(argument_map, print_usage_function):
        filename = opth.validate_required('input', argument_map, print_usage_function)
        output_filename = opth.validate_required('output', argument_map, print_usage_function)
        percent = opth.validate_required_float('percentage', argument_map, print_usage_function)

        seq_type = 'fasta'
        if file_tools.looks_fastq(filename):
            seq_type = 'fastq'

        with open(output_filename, 'w') as handle:
            def probability_write(handle, seq_rec):
                if random() <= percent:
                    SeqIO.write(seq_rec, handle, seq_type)
            write_function = functools.partial(probability_write, handle)
            apply_function_to_sequences(filename, write_function, seq_type)


# !---------- 'edit' command ----------!
class CommandEdit(object):
    @staticmethod
    def get_description():
        return 'Edit the sequence id names and write them to a new file.'

    @staticmethod
    def get_option_data():
        options = opth.default_option_map_input_output()
        options['input']['description'] = 'An input sequence file'
        options['output']['description'] = 'An output sequence file'
        options['map'] = {'order': 3,
                          'short': 'm',
                          'long': 'map',
                          'input_name': '<id_map>',
                          'description': "A tab delimited map file of the form 'old_id [tab] new_id.'",
                          'optional': True}
        return options

    @staticmethod
    def get_command_data(order=0):
        command_map = {'description': CommandEdit.get_description(),
                       'options': CommandEdit.get_option_data(),
                       'order': order}
        return command_map

    @staticmethod
    def run_program(argument_map, print_usage_function):
        input_filename = opth.validate_required('input', argument_map, print_usage_function)
        output_filename = opth.validate_required('output', argument_map, print_usage_function)

        seq_type = 'fasta'
        if file_tools.looks_fastq(input_filename):
            seq_type = 'fastq'

        map_mode = opth.has_option('map', argument_map)
        if map_mode:
            id_map = file_tools.parse_simple_map(argument_map['map'])

        if map_mode:
            def rename_and_write(output_handle, input_seq_type, seq_rec):
                current_id = seq_rec.id
                if current_id in id_map:
                    seq_rec.id = id_map[current_id]
                    SeqIO.write(seq_rec, output_handle, input_seq_type)
            apply_function_and_output(input_filename, output_filename, rename_and_write, seq_type)
        else:
            def edit_seq(output_handle, input_seq_type, seaq_rec):
                current_id = seaq_rec.id
                new_name = current_id
                response = input("##> Current sequence: %s. Would you like to rename it? (y/n) " % current_id)
                if response in ('y', 'Y', 'yes', 'YES'):
                    continue_edit = True
                    while continue_edit:
                        new_name = input("New id for %s: " % current_id)
                        finish_response = input("Confirm rename sequence %s to %s? (y/n) "
                                                    % (current_id, new_name))
                        if finish_response in ('y', 'Y', 'yes', 'YES'):
                            continue_edit = False
                seaq_rec.id = new_name
                SeqIO.write(seaq_rec, output_handle, input_seq_type)
            apply_function_and_output(input_filename, output_filename, edit_seq, seq_type)


# !---------- 'split' command ----------!
class CommandSplit(object):
    @staticmethod
    def get_description():
        return 'Split a genome file into separate sequence files.'

    @staticmethod
    def get_option_data():
        options = opth.default_option_map_input_output()
        options['input']['description'] = 'An input sequence file'
        options['output']['description'] = 'A directory where each sequence file will be placed.'
        options['output']['input_name'] = '<output_directory>'
        return options

    @staticmethod
    def get_command_data(order=0):
        command_map = {'description': CommandSplit.get_description(),
                       'options': CommandSplit.get_option_data(),
                       'order': order}
        return command_map

    @staticmethod
    def run_program(argument_map, print_usage_function):
        input_filename = opth.validate_required('input', argument_map, print_usage_function)
        output_directory = opth.validate_required('output', argument_map, print_usage_function)

        seq_type = 'fasta'
        if file_tools.looks_fastq(input_filename):
            seq_type = 'fastq'
            ext = '.fq'
        else:
            ext = '.fa'

        def write_seq_file(rec):
            seq_id = rec.id
            output_filename = file_tools.join(output_directory, seq_id + ext)
            with open(output_filename, 'w') as output_handle:
                SeqIO.write(rec, output_handle, seq_type)
        apply_function_to_sequences(input_filename, write_seq_file, seq_type)


# !---------- `prefix` command ----------!
class CommandPrefix(object):
    @staticmethod
    def get_description():
        return 'Add a prefix to every sequence id.'

    @staticmethod
    def get_option_data():
        options = opth.default_option_map_input_output()
        options['input']['description'] = 'An input sequence file'
        options['output']['description'] = 'An output sequence file'
        options['prefix'] = {'order': 3,
                             'short': 'p',
                             'long': 'prefix',
                             'input_name': '<prefix_text>',
                             'description': "The text to be used as a prefix and added to each record id",
                             'optional': False}
        return options

    @staticmethod
    def get_command_data(order=0):
        command_map = {'description': CommandPrefix.get_description(),
                       'options': CommandPrefix.get_option_data(),
                       'order': order}
        return command_map

    @staticmethod
    def run_program(argument_map, print_usage_function):
        input_filename = opth.validate_required('input', argument_map, print_usage_function)
        output_filename = opth.validate_required('output', argument_map, print_usage_function)
        prefix = opth.validate_required('prefix', argument_map, print_usage_function)

        seq_type = 'fasta'
        if file_tools.looks_fastq(input_filename):
            seq_type = 'fastq'

        def write_updated_seq_file(output_handle, input_seq_type, seq_rec):
            seq_id = seq_rec.id
            new_id = '_'.join([prefix, seq_id])
            seq_rec.id = new_id
            SeqIO.write(seq_rec, output_handle, input_seq_type)
        apply_function_and_output(input_filename, output_filename, write_updated_seq_file, seq_type)


# !---------- 'partition' command ----------!
class CommandPartition(object):
    @staticmethod
    def get_description():
        return 'Partition the sequences into n smaller files by a random process.'

    @staticmethod
    def get_option_data():
        options = opth.default_option_map_input_output()
        options['input']['description'] = 'An input sequence file'
        options['output']['description'] = 'the prefix to append to output files'
        options['output']['input_name'] = '<output_prefix>'
        options['number'] = {'order': 3,
                             'short': 'n',
                             'long': 'number',
                             'input_name': '<number_of_partitions>',
                             'description': "The number of partitions to be created.",
                             'optional': False}
        return options

    @staticmethod
    def get_command_data(order=0):
        command_map = {'description': CommandPartition.get_description(),
                       'options': CommandPartition.get_option_data(),
                       'order': order}
        return command_map

    @staticmethod
    def run_program(argument_map, print_usage_function):
        input_filename = opth.validate_required('input', argument_map, print_usage_function)
        output_prefix = opth.validate_required('output', argument_map, print_usage_function)
        number = opth.validate_required_int('number', argument_map, print_usage_function)

        seq_type = 'fasta'
        if file_tools.looks_fastq(input_filename):
            seq_type = 'fastq'

        files = []
        for i in range(number):
            output_file = open(output_prefix + '_' + str(i) + '.' + seq_type, 'w')
            files.append(output_file)

        def write_seq_to_random_file(seq):
            index = randrange(len(files))
            SeqIO.write(seq, files[index], seq_type)

        apply_function_to_sequences(input_filename, write_seq_to_random_file, seq_type)

        for output_file in files:
            output_file.close()


# !---------- base-comp command ----------!
class CommandBaseComposition(object):
    @staticmethod
    def get_description():
        return 'Calculate the base distribution composition for each sequence.'

    @staticmethod
    def get_option_data():
        options = opth.default_option_map_input()
        options['input']['description'] = 'An input fasta file'
        options['input']['input_name'] = '<sequence_file>'
        return options

    @staticmethod
    def get_command_data(order=0):
        command_map = {'description': CommandBaseComposition.get_description(),
                       'options': CommandBaseComposition.get_option_data(),
                       'order': order}
        return command_map

    @staticmethod
    def run_program(argument_map, print_usage_function):
        opth.validate_required('input', argument_map, print_usage_function)
        filename = argument_map['input']

        seq_type = 'fasta'
        if file_tools.looks_fastq(filename):
            seq_type = 'fastq'

        bases = ['A', 'C', 'G', 'T']
        header = ['sequence']
        header.extend(bases)
        print('\t'.join(header))

        def print_distribution_func(seq_rec):
            base_map = defaultdict(int)
            for base in str(seq_rec.seq).upper():
                base_map[base] += 1
            for base in str(seq_rec.seq.reverse_complement()).upper():
                base_map[base] += 1

            if len(base_map) > 5:
                print_usage_function('ERROR: Amino Acid sequences not yet supported for base-comp.')
                opth.sys_exit(1)

            output = [seq_rec.id]
            use_length = len(seq_rec.seq)
            if 'N' in base_map:
                use_length -= (base_map['N']/2)
            for base in bases:
                percent = base_map[base]*1.0/(use_length * 2)
                output.append(str(round(percent, 4)))
            print('\t'.join(output))
        apply_function_to_sequences(filename, print_distribution_func, seq_type)


# !---------- 'bed-fasta' command ----------!
class CommandBedFasta(object):
    @staticmethod
    def get_description():
        return 'Get the fasta sequences based on bed regions.'

    @staticmethod
    def get_option_data():
        options = opth.default_option_map_input_output()
        options['input']['description'] = 'An input sequence file'
        options['output']['description'] = 'The output fasta file with the extracted regions.'
        options['output']['input_name'] = '<output_file>'
        options['bed'] = {'order': 3,
                          'short': 'b',
                          'long': 'bed',
                          'input_name': '<bed_file>',
                          'description': "A 3 column bed file representing regions to extract.",
                          'optional': False}
        options['flank'] = {'order': 3,
                            'short': 'f',
                            'long': 'flank',
                            'input_name': '<flank_size>',
                            'description': "The size of a flank added to expand the bed region.",
                            'optional': True}
        options['prefix'] = {'order': 3,
                             'short': 'p',
                             'long': 'prefix',
                             'input_name': '<prefix_text>',
                             'description': "The text to be used as a prefix and added to each record id.",
                             'optional': True}
        return options

    @staticmethod
    def get_command_data(order=0):
        command_map = {'description': CommandBedFasta.get_description(),
                       'options': CommandBedFasta.get_option_data(),
                       'order': order}
        return command_map

    @staticmethod
    def run_program(argument_map, print_usage_function):
        filename = opth.validate_required('input', argument_map, print_usage_function)
        output_filename = opth.validate_required('output', argument_map, print_usage_function)
        bed_filename = opth.validate_required('bed', argument_map, print_usage_function)

        flank = opth.with_default_int('flank', argument_map, 0)
        prefix = opth.with_default('prefix', argument_map, '')

        if file_tools.looks_fastq(filename):
            print_usage_function("ERROR: bed-fasta requires that the input is a fasta file")
            opth.sys_exit(1)

        sequence_db = {}
        with open(filename) as f:
            for seq in SeqIO.parse(f, 'fasta'):
                sequence_db[seq.id] = seq

        fout = open(output_filename, 'w')
        with open(bed_filename) as f:
            for line in f:
                line = line.strip()
                rec = line.split('\t')
                rec_id = rec[0]
                start = int(rec[1])
                end = int(rec[2])

                if rec_id in sequence_db:
                    seq = sequence_db[rec_id]
                    new_start = start - flank
                    if new_start < 0:
                        new_start = 0
                    new_end = end + flank
                    if new_end > len(seq):
                        new_end = len(seq)

                    if len(prefix) > 0:
                        new_prefix = prefix + '|'
                    else:
                        new_prefix = ''

                    new_id = new_prefix + rec_id + ':' + str(new_start) + '-' + str(new_end)
                    new_rec = SeqRecord.SeqRecord(Seq.Seq(str(seq.seq[new_start:new_end])), id=new_id, description='')
                    SeqIO.write(new_rec, fout, 'fasta')
        fout.close()


# !---------- 'to-fasta' command ----------!
class CommandToFasta(object):
    @staticmethod
    def get_description():
        return 'Convert fastq files to fasta files.'

    @staticmethod
    def get_option_data():
        options = opth.default_option_map_input_output()
        options['input']['description'] = 'An input fastq file'
        options['input']['input_name'] = '<input_file>'
        options['output']['description'] = 'The output fasta file.'
        options['output']['input_name'] = '<output_file>'
        return options

    @staticmethod
    def get_command_data(order=0):
        command_map = {'description': CommandToFasta.get_description(),
                       'options': CommandToFasta.get_option_data(),
                       'order': order}
        return command_map

    @staticmethod
    def run_program(argument_map, print_usage_function):
        filename = opth.validate_required('input', argument_map, print_usage_function)
        output_filename = opth.validate_required('output', argument_map, print_usage_function)

        def write_sequence(output_handle, seq_type, seq_record):
            SeqIO.write(seq_record, output_handle, 'fasta')
        apply_function_and_output(filename, output_filename, write_sequence, 'fastq')


# !---------- `number` command ----------!
class CommandNumber(object):
    @staticmethod
    def get_description():
        return 'Add an incrementing number of every sequence ID.'

    @staticmethod
    def get_option_data():
        options = opth.default_option_map_input_output()
        options['input']['description'] = 'An input sequence file'
        options['output']['description'] = 'An output sequence file'
        options['prefix'] = {'order': 3,
                             'short': 'p',
                             'long': 'prefix',
                             'input_name': '<prefix_text>',
                             'description': "The text to be used as a prefix and added to each record id",
                             'optional': True}
        return options

    @staticmethod
    def get_command_data(order=0):
        command_map = {'description': CommandNumber.get_description(),
                       'options': CommandNumber.get_option_data(),
                       'order': order}
        return command_map

    @staticmethod
    def run_program(argument_map, print_usage_function):
        input_filename = opth.validate_required('input', argument_map, print_usage_function)
        output_filename = opth.validate_required('output', argument_map, print_usage_function)
        prefix = opth.with_default('prefix', argument_map, '')

        seq_type = 'fasta'
        if file_tools.looks_fastq(input_filename):
            seq_type = 'fastq'

        seq_count = 0
        output_handle = open(output_filename, 'w')
        with open(input_filename) as handle:
            for seq_rec in SeqIO.parse(handle, seq_type):
                seq_id = seq_rec.id
                new_id = '_'.join([prefix, seq_id, str(seq_count)])
                seq_count += 1
                seq_rec.id = new_id
                SeqIO.write(seq_rec, output_handle, seq_type)
        output_handle.close()


# !---------- 'sliding' command ----------!
class CommandSliding(object):
    @staticmethod
    def get_description():
        return 'Divide a sequence into subsequences using a sliding window method.'

    @staticmethod
    def get_option_data():
        options = opth.default_option_map_input_output()
        options['input']['description'] = 'An input sequence file'
        options['output']['description'] = 'An output sequence file'
        options['window'] = {'order': 3,
                             'short': 'w',
                             'long': 'window',
                             'input_name': '<window_size>',
                             'description': "The size of the window in bp.",
                             'optional': False}
        options['step'] = {'order': 4,
                           'short': 's',
                           'long': 'step',
                           'input_name': '<step_size>',
                           'description': "The size of the step between window sequences.",
                           'optional': False}
        return options

    @staticmethod
    def get_command_data(order=0):
        command_map = {'description': CommandSliding.get_description(),
                       'options': CommandSliding.get_option_data(),
                       'order': order}
        return command_map

    @staticmethod
    def run_program(argument_map, print_usage_function):
        input_filename = opth.validate_required('input', argument_map, print_usage_function)
        output_filename = opth.validate_required('output', argument_map, print_usage_function)
        window_size = opth.validate_required_int('window', argument_map, print_usage_function)
        step_size = opth.validate_required_int('step', argument_map, print_usage_function)

        seq_type = 'fasta'
        if file_tools.looks_fastq(input_filename):
            seq_type = 'fastq'

        def gen_window_seqs(output_handle, input_seq_type, seq_rec):
            seq_id = seq_rec.id
            seq_len = len(seq_rec)
            start = 0
            end = start + window_size
            while start < seq_len - window_size:
                new_seq = seq_rec[start:end]
                new_id = '_'.join([seq_id,str(start),str(end)])
                new_seq.id = new_id
                SeqIO.write(seq_rec, output_handle, input_seq_type)

                # slide window by step
                start = start + step_size
                end = start + window_size

        apply_function_and_output(input_filename, output_filename, gen_window_seqs, seq_type)
