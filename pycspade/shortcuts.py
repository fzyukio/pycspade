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
    def __init__(self, name, xname, yname, noccurs):
        self.items = []
        self.name = name
        self.xname = xname
        self.yname =yname
        self.noccurs = noccurs
        self.confidence = None
        self.lift = None

    def add_item(self, item):
        self.items.append(item)

    def __repr__(self):
        return '{} - [{}]'.format('->'.join(list(map(str, self.items))), self.noccurs)


def decode_results(result):
    occurrences = {}
    lifts = {}
    confidences = {}
    nseqs = result['nsequences']

    mined = result['mined']
    lines = mined.strip().decode('latin-1').split('\n')
    lines.sort()
    sequences = []
    for line in lines:
        if '0' <= line[0] <= '9':
            _sequence, stats = line.split(' -- ')
            _items = _sequence.split(' -> ')
            noccurs = int(stats[:stats.index(' ')])
            # lift = None

            if len(_items) > 1:
                x_item = ' -> '.join(_items[:-1])
                y_item = _items[-1]
            else:
                x_item = None
                y_item = None

            occurrences[_sequence] = noccurs
            sequence = Sequence(_sequence, x_item, y_item, noccurs)

            for _item in _items:
                _elements = list(map(int, _item.split(' ')))
                item = Item(_elements)
                sequence.add_item(item)
            sequences.append(sequence)

    # Second pass
    for sequence in sequences:
        if sequence.lift is not None:
            continue
        x_item = sequence.xname
        y_item = sequence.yname

        x_noccurs = occurrences.get(x_item, None)
        y_noccurs = occurrences.get(y_item, None)

        if x_noccurs is not None:
            sequence.confidence = sequence.noccurs / x_noccurs
            confidences[sequence.name] = sequence.confidence

            if y_noccurs is not None:
                sequence.lift = sequence.noccurs * nseqs / (x_noccurs * y_noccurs)
                lifts[sequence.name] = sequence.lift

    result['mined_objects'] = sequences


def cspade(filename=None, data=None, support=3, maxsize=None, maxlen=None, mingap=None, maxgap=None):
    """
    Shortcut to call cspade
    :param filename: path to the ascii file, must be given if data is None
    :param data: raw data as list of transactions, must be given if filename is None
    :param support: is interpreted as the threshold of mimimum normalised support if within [0, 1]:
                         if > 1: interpreted as the threshold of absolute support (e.g. 50 over 100 transactions)
    :param maxsize: an integer value specifying the maximum number of items of a sequence (default=100)
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
        return retval
    finally:
        if data:
            os.remove(filename)
