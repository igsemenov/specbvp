# -*- coding: utf-8 -*-
"""Documentation maker.
"""
import runpy
import docspyer

SRCPATH = 'docs/sources'
DOCPATH = 'docs/build'

config = {
    'swaplinks': True
}

runpy.run_path('_docmods.py')

docspyer.builddocs(
    SRCPATH, DOCPATH, **config
)
