import os
import uuid

from libcpp cimport bool
from libcpp.string cimport string as c_string
from libcpp.vector cimport vector

cdef extern from "../csrc/sequence.cc":
    vector[c_string] mine(char *dbname, double min_support_one, int min_support_all, int use_ascending, bool use_class,
                int _num_partitions, bool ext_l2_pass, bool use_hash, int min_gap, int avaimem_mb, bool _outputfreq,
                bool recursive, int pruning_type, int max_gap, bool use_window, int max_seq_len, int max_iset_len
                )  except +


cdef extern from "../csrc/makebin.cc":
    void convert_bin(c_string ifname, c_string ofname)  except +

cdef extern from "../csrc/getconf.cc":
    c_string create_conf(c_string datafile_name, bool assoc) except +

cdef extern from "../csrc/exttpose.cc":
    vector[c_string] _exttpose(c_string dbname, int num_partitions, double min_support, bool twoseq, bool use_diff, bool do_l2,
                  bool do_invert, bool use_newformat, int maxmem, bool no_minus_off) except +


def cpp_cspade(filename, decode=True, num_partitions = 1, min_support_one = 0., min_support_all = -1, twoseq = True,
                 use_diff = False, do_l2 = False, do_invert = True, use_newformat = True, maxmem = 128,
                 no_minus_off = True, use_ascending = -2, use_class = False, ext_l2_pass = True,
                 use_hash = False, min_gap = 1, recursive = True, pruning_type = 0, max_gap = 1, use_window = False,
                 max_seq_len = 100, max_iset_len = 100):

    avaimem_mb = 128
    outputfreq = True


    dbname = uuid.uuid4().hex
    binfile = bytes('/tmp/{}.data'.format(dbname), encoding='latin-1')
    dbname = bytes('/tmp/{}'.format(dbname), encoding='latin-1')

    convert_bin(bytes(filename, encoding='latin-1'), binfile)
    conffile = create_conf(dbname, False)

    tposefiles = _exttpose(dbname, num_partitions, min_support_one, twoseq, use_diff, do_l2, do_invert, use_newformat,
                           maxmem, no_minus_off)
    seq, log = mine(dbname, min_support_one, min_support_all, use_ascending, use_class, num_partitions, ext_l2_pass,
               use_hash, min_gap, avaimem_mb, outputfreq, recursive, pruning_type, max_gap, use_window, max_seq_len,
               max_iset_len)

    tmp_files = tposefiles + [conffile, binfile]

    for tmp_file in tmp_files:
        if os.path.isfile(tmp_file):
            os.remove(tmp_file)

    if decode:
        seq = seq.decode('latin-1').split('\n')
        log = log.decode('latin-1').split('\n')

    return seq, log
