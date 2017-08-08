# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""C libraries for miRmap."""

from . import if_clib_rna
from . import if_clib_spatt
from . import if_clib_phast

class LibraryLink(object):
    def __init__(self, library_path=None):
        """:param library_path: Path to the C dynamic libraries.
           :type library_path: str"""
        self.libs = {}
        if library_path is not None:
            self.init_library_link('rna', library_path=library_path)
            self.init_library_link('spatt', library_path=library_path)
            self.init_library_link('phast', library_path=library_path)

    def init_library_link(self, library_short_name, library_path=None, library_name=None):
        if library_short_name == 'rna':
            self.libs[library_short_name] = if_clib_rna.RNAvienna(library_path, library_name)
        if library_short_name == 'spatt':
            self.libs[library_short_name] = if_clib_spatt.Spatt(library_path, library_name)
        if library_short_name == 'phast':
            self.libs[library_short_name] = if_clib_phast.Phast(library_path, library_name)

    def get_library_link(self, library_short_name):
        try:
            return self.libs[library_short_name]
        except KeyError:
            self.init_library_link(library_short_name)
            return self.libs[library_short_name]
