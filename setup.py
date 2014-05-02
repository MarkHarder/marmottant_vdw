from distutils.core import setup
import io
import os
import sys

import marmottant_vdw.marmottant_vdw

here = os.path.abspath(os.path.dirname(__file__))

def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(filename, encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.txt')

class PyTest():
    def finalize_options(self):
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        errcode = pytest.main(self.test_args)
        sys.exit(errcode)

setup(
    name='marmottant_vdw',
    version=marmottant_vdw.marmottant_vdw.__version__,
    url='http://github.com/markharder/marmottant_vdw/',
    author='Mark Harder',
    cmdclass={'test': PyTest},
    description='Marmottant Van Der Waal bubble simulation',
    long_description=long_description,
    packages=['marmottant_vdw'],
    platforms='any',
    classifiers = [
        'Programming Language :: Python',
        'Development Status :: 4 - Beta',
        'Natural Language :: English',
        'Intended Audience :: Academic',
        'Operating System :: OS Independent',
        ],
)
