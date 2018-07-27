import os
import uuid

from libcpp cimport bool
from libcpp.string cimport string as c_string

cdef extern from "../csrc/utils.h":
    cdef enum:
        Pruning_No = 0
        Pruning_L2 = 1
        Pruning_Zero = 2
        Pruning_Follow = 4

cdef extern from "../csrc/utils.h":
    cdef struct result_t:
        int nsequences;
        bool success;
        c_string error;
        c_string mined;
        c_string logger;
        c_string summary;
        c_string memlog;


cdef extern from "../csrc/utils.h" namespace "sequence":
    cdef struct arg_t:
        c_string name
        c_string binf
        c_string dataf
        c_string idxf
        c_string conf
        c_string it2f
        c_string seqf
        c_string classf

        int num_partitions
        double min_support
        int use_ascending
        bool use_class
        bool do_l2
        bool use_hash
        int min_gap
        double maxmem
        bool recursive
        int pruning_type
        int max_gap
        int max_seq_len
        int max_iset_len

        bool twoseq
        bool use_diff
        bool use_newformat
        bool no_minus_off


cdef extern from "../csrc/sequence.cc":
    result_t cspade(c_string asciifile, arg_t& _args)


cdef extern from "../csrc/sequence.cc" namespace "sequence":
    arg_t populate_names(arg_t &_args)

def cpp_cspade(filename, support=3, maxsize=None, maxlen=None, mingap=None, maxgap=None, decode=True):
    """
    Call C++'s cspade()
    :param filename: full path to the input file (ascii)
    :param support: is interpreted as the threshold of mimimum normalised support if within [0, 1]:
                         if > 1: interpreted as the threshold of absolute support (e.g. 50 over 100 transactions)
    :param maxsize: an integer value specifying the maximum number of items of an element of a sequence (default=100)
    :param maxlen: an integer value specifying the maximum number of elements of a sequence (default=100)
    :param mingap: an integer value specifying the minimum time difference between consecutive elements of a sequence
    :param maxgap: an integer value specifying the minimum time difference between consecutive elements of a sequence
    :param decode: if True, the return strings will be decoded and line-separated, otherwise raw C++ strings
                   (python bytes) are returned
    :return: (result, logger, summary, memlog). where:
             -result: the mined sequences
             -logger: general logging
             -summary: equivalent to the content of summary.out
             -memlog: logging of memory usage
    """
    assert (support > 0 and (support < 1 or (support >= 1 and float(support).is_integer()))), \
           'support must be a floating point in range [0-1] (percentage) or an int >= 1 (absolute)'
    avaimem_mb = 128

    dbname = uuid.uuid4().hex
    dbname = bytes('/tmp/{}'.format(dbname), encoding='latin-1')
    filename = bytes(filename, encoding='latin-1')

    cdef arg_t args
    args.num_partitions = 1
    args.min_support = support
    args.use_ascending = -2
    args.use_class = False
    args.do_l2 = False
    args.use_hash = False
    args.min_gap = 1
    args.max_gap = 2147483647
    args.maxmem = 128
    args.recursive = False
    args.pruning_type = Pruning_No

    args.max_seq_len = 100
    args.max_iset_len = 100

    args.twoseq = False
    args.use_diff = False
    args.use_newformat = True
    args.no_minus_off = True

    args.name = dbname

    if maxlen is not None:
        args.max_seq_len = maxlen
    if maxsize is not None:
        args.max_iset_len = maxsize
    if mingap is not None:
        assert mingap > 0, 'mingap cannot be 0 - that would mean two transactions happen at the same time'
        args.min_gap = mingap
    if maxgap is not None:
        assert maxgap > 0, 'maxgap cannot be 0'
        args.max_gap = maxgap
        if args.max_gap < args.min_gap:
            args.min_gap = args.max_gap

    args = populate_names(args)

    cdef result = cspade(filename, args)

    tmp_files = [args.binf, args.conf, args.idxf, args.it2f, args.dataf, args.seqf, args.classf]
    for tmp_file in tmp_files:
        if os.path.isfile(tmp_file):
            os.remove(tmp_file)

    if result['success']:
        if decode:
            result['mined'] = result['mined'].decode('latin-1').split('\n')
            result['logger'] = result['logger'].decode('latin-1').split('\n')
            result['summary'] = result['summary'].decode('latin-1').split('\n')
            result['memlog'] = result['memlog'].decode('latin-1').split('\n')
        return result
    else:
        raise RuntimeError(result['error'].decode('latin-1'))