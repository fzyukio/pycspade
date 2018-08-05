from setuptools import setup, Extension
from codecs import open

try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

if use_cython:
    sourcefiles = ['pycspade/cspade.pyx']
else:
    sourcefiles = ['pycspade/cspade.cpp']

other_files = ['csrc/{}'.format(x) for x in [
    'Itemset.cc', 'Array.cc', 'ArrayT.cc', 'Eqclass.cc', 'Lists.cc', 'extl2.cc', 'partition.cc', 'maxgap.cc',
    'calcdb.cc', 'makebin.cc', 'getconf.cc', 'exttpose.cc', 'utils.cc'
]]

sourcefiles += other_files

ext_modules = [
    Extension('pycspade.cspade',
              sourcefiles,
              include_dirs=['csrc/'],
              language='c++',
              extra_compile_args=[
                  '-std=c++11',
                  '-Wno-sign-compare',
                  '-Wno-incompatible-pointer-types',
                  '-Wno-unused-variable',
                  '-Wno-absolute-value',
                  '-Wno-visibility',
                  '-Wno-#warnings']
              )
]

setup(
    name='pycspade',
    cmdclass={'build_ext': build_ext},
    ext_modules=ext_modules,
    license='MIT',
    packages=['pycspade'],
    version='0.3.2',
    author=['Mohammed J. Zaki', 'Yukio Fukuzawa'],
    description='C-SPADE Python Implementation',
    long_description=open('README.md').read(),
    url='https://github.com/fzyukio/pycspade',
    keywords=['cspade', 'c-spade', 'sequence mining'],
    install_requires=['Cython'],

)
