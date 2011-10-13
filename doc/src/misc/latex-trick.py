#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Modify sphinx-produced latex directory for Japanese documents

 $Id$
"""

import os
import sys
import re
import string
import shutil
import codecs
import tempfile

THISDIR   = os.path.dirname(os.path.abspath(__file__))
SRCDIR    = string.join(THISDIR.split('/')[:-1], '/')
sys.path.append(SRCDIR)
import conf

LATEXSRC  = conf.latex_documents[0][1]
MAKEFILE  = os.path.join(THISDIR, 'Makefile')
SPHINXSTY = os.path.join(THISDIR, 'sphinx.sty')
REPLACE_DICT = {
    r'^(\\usepackage{babel})' : r'%%%\1',
}

def replace(filename, dic):
    if not os.path.exists(filename):
        raise ValueError, '%s does not exist !' % filename
    fp = codecs.open(filename, 'r', 'utf-8')
    try:
        tempname = tempfile.mkstemp(suffix='.tex', prefix='temp-')[1]
        tf = open(tempname, 'w')
        l = fp.readline()
        while l != '':
            for key, item in dic.items():
                r = re.sub(key, item, l)
                l = r
            # encoding needed for legacy platex
            tf.write(l.encode('euc_jp'))
            l = fp.readline()
        tf.close()
        shutil.copy(tempname, filename)
    finally:
        os.remove(tempname)

if __name__ == '__main__':
    latex_dir = sys.argv[1]
    if not os.path.exists(latex_dir):
        raise ValueError, '%s directory does not exists' % latex_dir
    # copy Makefile
    shutil.copy(MAKEFILE, os.path.join(latex_dir, 'Makefile'))
    # copy sphinx.sty
    shutil.copy(SPHINXSTY, os.path.join(latex_dir, 'sphinx.sty'))
    # modify tex source
    replace(os.path.join(latex_dir, LATEXSRC), REPLACE_DICT)
