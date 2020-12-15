"""
Some simple functions to work with files and file names
"""
import os


def curdir():
    return os.path.curdir


def abs_curdir():
    return os.path.abspath(os.path.curdir)


def basename(filename):
    return os.path.basename(filename)


def abspath(filename):
    return os.path.abspath(filename)


def parse_lines_to_set(filename):
    line_set = set()
    with open(filename) as f:
        for line in f:
            line_set.add(line.strip())
    return line_set


def parse_lines_to_list(filename):
    line_list = []
    with open(filename) as f:
        for line in f:
            line_list.append(line.strip())
    return line_list


def parse_simple_map(filename, separator='\t'):
    line_map = {}
    with open(filename) as f:
        for line in f:
            record = line.split(separator)
            line_map[record[0]] = record[1]
    return line_map


def fastq_base_name(string):
    string = os.path.basename(string)
    if string.endswith('.gz'):
        tmp_string = string[:-3]
    else:
        tmp_string = string

    if tmp_string.endswith('.fq'):
        tmp_string = tmp_string[:-3]
    elif tmp_string.find('.fastq'):
        tmp_string = tmp_string[:-6]
    return tmp_string


def directory_exists(in_dir):
    return os.path.isdir(in_dir)


def looks_like_maker(filename_suffix):
    def looks_like(filename):
        index = filename.find(filename_suffix)
        return index >= 0 and (len(filename) - index) == len(filename_suffix)
    return looks_like


def looks_zipped(filename):
    looks_like = looks_like_maker('.gz')
    return looks_like(filename)


def looks_fastq(filename):
    looks_like_fastq = looks_like_maker('.fastq')
    looks_like_fq = looks_like_maker('.fq')
    return looks_like_fastq(filename) or looks_like_fq(filename)


def trim_maker(filename_siffix):
    def trim_suffix(filename):
        index = filename.find(filename_siffix)
        if index > 0:
            return filename[:index]
        else:
            return filename
    return trim_suffix


def shared_prefix(input_filename1, input_filename2):
    filename_1 = basename(input_filename1)
    filename_2 = basename(input_filename2)

    min_len = min(len(filename_1), len(filename_2))
    i = 0
    while i < min_len and (filename_1[i] == filename_2[i]):
        i += 1
    return filename_1[:i].strip('._')


def join(path, *paths):
    return os.path.join(path, *paths)


def trim_pair_number(read_id):
    trim_1 = trim_maker('/1')
    trim_2 = trim_maker('/2')
    return trim_1(trim_2(read_id))

