import sys
print (sys.path)
sys.path.insert(0, "pycspade")

try:
    from pycspade import cspade, cpp_cspade
except ImportError:
    print('No Import')

# seq, log, summary, memlog = cpp_cspade('tests/test1.ascii.data', 0.2, 100, decode=False)
# print(seq.decode('latin-1'))

result = cspade(filename='tests/test1.ascii.data', support=20, maxsize=5, maxlen=5)
nseqs = result['nsequences']

print('{0:>9s} {1:>9s} {2:>9s} {3:>9s} {4:>80s}'.format('Occurs', 'Support', 'Confid', 'Lift', 'Sequence'))
for mined_object in result['mined_objects']:
    conf = 'N/A'
    lift = 'N/A'
    if mined_object.confidence:
        conf = '{:0.7f}'.format(mined_object.confidence)
    if mined_object.lift:
        lift = '{:0.7f}'.format(mined_object.lift)

    print('{0:>9d} {1:>0.7f} {2:>9s} {3:>9s} {4:>80s} '.format(
        mined_object.noccurs, mined_object.noccurs / nseqs, conf, lift, '->'.join(list(map(str, mined_object.items)))))


# data = [
#     [1, 1, [1]],
#     [1, 2, [2]],
#     [1, 3, [3]],
#
#     [2, 1, [2]],
#     [2, 2, [3]],
#
#     [3, 1, [1]],
#     [3, 2, [3]],
#     [3, 3, [2]],
#     [3, 4, [3]],
# ]
#
# seq, log, summary, memlog = cspade(data=data, support=0.5, decode=False)
# print(seq.decode('latin-1'))