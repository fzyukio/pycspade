import sys
print (sys.path)
sys.path.insert(0, "pycspade")

try:
    from pycspade import cspade, cpp_cspade
except ImportError:
    print('No Import')

# seq, log, summary, memlog = cpp_cspade('tests/test1.ascii.data', 0.2, 100, decode=False)
# print(seq.decode('latin-1'))

result = cspade(filename='tests/test1.ascii.data', support=100, maxgap=1)
for mined_object in result['mined_objects']:
    print(mined_object)


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