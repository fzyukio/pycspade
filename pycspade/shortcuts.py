import os
import uuid

from .cspade import cpp_cspade


def data_to_rows(data):
    rows = ['{} {} {} {}'.format(sid, eid, len(els), ' '.join(list(map(str, els)))) for sid, eid, els in data]
    return rows


class Item:
    def __init__(self, elements):
        self.elements = elements

    def __repr__(self):
        return '({})'.format(' '.join(list(map(str, self.elements))))


class Sequence:
    def __init__(self, support):
        self.items = []
        self.support = support

    def add_item(self, item):
        self.items.append(item)

    def __repr__(self):
        return '{} - [{}]'.format('->'.join(list(map(str, self.items))), self.support)


def decode_results(result):
    mined = result['mined']
    nseqs = result['nsequences']
    lines = mined.strip().decode('latin-1').split('\n')
    sequences = []
    for line in lines:
        if '0' <= line[0] <= '9':
            _sequence, stats = line.split(' -- ')
            _items = _sequence.split(' -> ')
            _support = int(stats[:stats.index(' ')])

            sequence = Sequence(_support)
            for _item in _items:
                _elements = list(map(int, _item.split(' ')))
                item = Item(_elements)
                sequence.add_item(item)
            sequences.append(sequence)
    result['mined_objects'] = sequences


def cspade(filename=None, data=None, support=3, maxsize=None, maxlen=None, mingap=None, maxgap=None):
    """
    Shortcut to call cspade
    :param filename: path to the ascii file, must be given if data is None
    :param data: raw data as list of transactions, must be given if filename is None
    :param support: is interpreted as the threshold of mimimum normalised support if within [0, 1]:
                         if > 1: interpreted as the threshold of absolute support (e.g. 50 over 100 transactions)
    :param maxsize: an integer value specifying the maximum number of items of an element of a sequence (default=100)
    :param maxlen: an integer value specifying the maximum number of elements of a sequence (default=100)
    :param mingap: an integer value specifying the minimum time difference between consecutive elements of a sequence
    :param maxgap: an integer value specifying the maximum time difference between consecutive elements of a sequence
    :param decode: if True, the return strings will be decoded and line-separated, otherwise raw C++ strings
                   (python bytes) are returned
    :return: (result, logger, summary, memlog). where:
             -result: the mined sequences
             -logger: general logging
             -summary: equivalent to the content of summary.out
             -memlog: logging of memory usage
    """
    if filename is None and data is None:
        raise Exception('You must provide either filename or data')
    if filename is not None and data is not None:
        raise Exception('You must provide either filename or data')

    if data:
        rows = data_to_rows(data)
        hex = uuid.uuid4().hex
        filename = '/tmp/{}.ascii.data'.format(hex)
        with open(filename, 'w', encoding='latin-1') as f:
            for row in rows:
                f.write(row)
                f.write('\n')

    try:
        retval = cpp_cspade(filename, support, maxsize, maxlen, mingap, maxgap, decode=False)
        decode_results(retval)
    finally:
        if data:
            os.remove(filename)

    return retval
