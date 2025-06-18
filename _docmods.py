# -*- coding: utf-8 -*-
"""Documentation maker — modules.
"""
import docspyer
import specbvp

DOCPATH = 'docs/sources'

MODULES = [
    specbvp.polybases,
]

config = {
    'docsname': 'specbvp',
    'hostname': 'specbvp',
    'modrefs': True,
    'clsverbs': 2
}

docspyer.docmods(
    MODULES, DOCPATH, **config
)
