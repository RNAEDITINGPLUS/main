# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Ctypes interface classe with the `Spatt <http://www.mi.parisdescartes.fr/~nuel/spatt>`_ C library."""

import os

from ctypes import *

from . import utils

class Spatt(object):
    """Interface class for the Spatt library."""
    def __init__(self, library_path=None, library_name=None):
        if library_name is None:
            lib = 'libspatt2.so'
        else:
            lib = library_name
        if library_path is not None:
            lib = os.path.join(library_path, lib)
        self._library = cdll.LoadLibrary(lib)
        # Functions arguments and result types
        self._library.spatt_exact.argtypes = [c_char_p, c_char_p, POINTER(c_double), c_short, c_bool, c_long, c_ulong, c_char]
        self._library.spatt_exact.restype = c_double

    def transitions_l2c(self, transitions):
        """Take the transition matrix a list with element in C-order"""
        ctransitions = (c_double * len(transitions))()
        i = 0
        for v in transitions:
            ctransitions[i] = v
            i += 1
        return ctransitions

    def get_exact_prob(self, motif, nobs, length_seq, alphabet, transitions, markov_order, direction):
        ctransitions = self.transitions_l2c(utils.flatten(transitions))
        return self._library.spatt_exact(''.join(alphabet), motif, cast(ctransitions, POINTER(c_double)), markov_order, False, nobs, length_seq, direction)
