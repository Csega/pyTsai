#!/usr/bin/env python

from distutils.core import *

pytsai_ext = Extension(
        'pytsai', [
        'src/pytsai.c',
        'src/errors.c',
        'src/tsai/cal_eval.c',
        'src/tsai/cal_main.c',
        'src/tsai/cal_tran.c',
        'src/tsai/ecalmain.c',
        'src/minpack/dpmpar.c',
        'src/minpack/enorm.c',
        'src/minpack/fdjac2.c',
        'src/minpack/lmdif.c',
        'src/minpack/lmpar.c',
        'src/minpack/qrfac.c',
        'src/minpack/qrsolv.c',
        'src/matrix/matrix.c'
])
#extra_compile_args=['-O2', '-Wall', '-pedantic', '-std=c99',
#'-W', '-Wunreachable-code'])

setup(
        name = 'PyTsai',
        version = '1.0',
        description = 'Tsai camera calibration.',
        author = 'Jonathan Merritt, Gabor Cseh',
        author_email = 'j.merritt@pgrad.unimelb.edu.au, csega@mailbox.hu',
        keywords = "tsai automatic calibration python blender",
		url = "https://github.com/Csega/pyTsai",   # project home page, if any
        ext_modules = [ pytsai_ext ],
        py_modules = [ 'Tsai' ]
)
