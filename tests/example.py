from pycspade import cspade

seq, log = cspade('tests/test.ascii.data', decode=False, min_support_one=0., max_gap=2)

print(seq)