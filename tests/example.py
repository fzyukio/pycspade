from pycspade import cspade

seq = cspade('tests/test.ascii.data', min_support_one=0., max_gap=2)

print(seq.decode('latin-1'))