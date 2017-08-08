# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Ctypes interface classes with the `Vienna RNA <http://www.tbi.univie.ac.at/RNA>`_ C library."""

import math
import os

from ctypes import *

# Classes for structures

class DUPLEX(Structure):
    _fields_ = [('i', c_int), ('j', c_int), ('structure', c_char_p), ('energy', c_float)]

class PU_CONTRIB(Structure):
    _fields_ = [('H', POINTER(POINTER(c_double))), ('I', POINTER(POINTER(c_double))), ('M', POINTER(POINTER(c_double))), ('E', POINTER(POINTER(c_double))), ('length', c_int)]

class PLIST(Structure):
    _fields_ = [('i', c_int), ('j', c_int), ('p', c_float)]

# - part_func_co.h

class COFOLDF(Structure):
    _fields_ = [('F0AB', c_double), ('FAB', c_double), ('FcAB', c_double), ('FA', c_double), ('FB', c_double)]

class RNAvienna(object):
    """Interface class for RNA library from Vienna"""
    def __init__(self, library_path=None, library_name=None):
        if library_name is None:
            lib = 'libRNAvienna.so'
        else:
            lib = library_name
        if library_path is not None:
            lib = os.path.join(library_path, lib)
        self._library = cdll.LoadLibrary(lib)
        # Functions arguments and result types
        self._library.fold.argtypes = [c_char_p, c_char_p]
        self._library.fold.restype = c_float
        self._library.Lfold.argtypes = [c_char_p, c_char_p, c_int]
        self._library.Lfold.restype = c_float
        self._library.init_pf_fold.argtypes = [c_int]
        self._library.init_pf_fold.restype = None
        self._library.pf_fold.argtypes = [c_char_p, c_char_p]
        self._library.pf_fold.restype = c_float
        self._library.pfl_fold.argtypes = [c_char_p, c_int, c_int, c_float, POINTER(c_double)]
        self._library.pfl_fold.restype = POINTER(PLIST)
        self._library.free_pf_arrays.argtypes = []
        self._library.free_pf_arrays.restype = None
        self._library.cofold.argtypes = [c_char_p, c_char_p]
        self._library.cofold.restype = c_float
        self._library.duplexfold.argtypes = [c_char_p, c_char_p]
        self._library.duplexfold.restype = DUPLEX
        self._library.pf_unstru.argtypes = [c_char_p, c_int]
        self._library.pf_unstru.restype = POINTER(PU_CONTRIB)
        self._library.space.argtypes = [c_int]
        self._library.space.restype = c_void_p
        self._library.mean_bp_distance.argtypes = [c_int]
        self._library.mean_bp_distance.restype = c_double
        # - part_func_co.h
        self._library.co_pf_fold.argtypes = [c_char_p, c_char_p]
        self._library.co_pf_fold.restype = COFOLDF
        self._library.free_co_pf_arrays.argtypes = []
        self._library.free_co_pf_arrays.restype = None
        self._library.init_co_pf_fold.argtypes = [c_int]
        self._library.init_co_pf_fold.restype = None

    def lfold(self, seq, struc_buffer, maxdist):
        return self._library.Lfold(seq, struc_buffer, maxdist)

    def duplexfold(self, seq1, seq2):
        return self._library.duplexfold(seq1, seq2)

    def init_pf_fold(self, length):
        self._library.init_pf_fold(length)

    def pf_fold(self, seq, struc_buffer):
        return self._library.pf_fold(seq, struc_buffer)

    def pfl_fold(self, seq, size_window, size_max_interact, size_free, ener_cutoff):
        pup = cast(self.space((len(seq)+1)*sizeof(c_double)), POINTER(c_double))
        if size_free > 0:
            pup[0] = size_free
        return self._library.pfl_fold(seq, size_window, size_max_interact, ener_cutoff, pup)

    def pfl_fold_with_pup(self, seq, size_window, size_max_interact, size_free, ener_cutoff):
        pup = cast(self.space((len(seq)+1)*sizeof(c_double)), POINTER(c_double))
        if size_free > 0:
            pup[0] = size_free
        pl = self._library.pfl_fold(seq, size_window, size_max_interact, ener_cutoff, pup)
        return pup, pl

    def free_pf_arrays(self):
        self._library.free_pf_arrays()

    def free_co_pf_arrays(self):
        self._library.free_co_pf_arrays()

    def pf_unstru(self, seq, length):
        return self._library.pf_unstru(seq, length)

    def init_co_pf_fold(self, length):
        self._library.init_co_pf_fold(length)

    def co_pf_fold(self, seq, struc_buffer):
        return self._library.co_pf_fold(seq, struc_buffer)

    def set_temperature(self, temperature):
        c_double.in_dll(self._library, 'temperature').value = temperature

    def set_dangles(self, dangles):
        c_int.in_dll(self._library, 'dangles').value = dangles

    def set_nolonelypairs(self, nolonelypairs):
        c_int.in_dll(self._library, 'noLonelyPairs').value = nolonelypairs

    def set_pf_scale(self, pf_scale):
        c_double.in_dll(self._library, 'pf_scale').value = pf_scale

    def set_fold_constrained(self, fold_constrained):
        c_int.in_dll(self._library, 'fold_constrained').value = fold_constrained

    def set_cut_point(self, cut_point):
        c_int.in_dll(self._library, 'cut_point').value = cut_point

    def get_pf_scale(self):
        return c_double.in_dll(self._library, 'pf_scale').value

    def get_iindx(self):
        return POINTER(c_int).in_dll(self._library, 'iindx')

    def get_pr(self):
        return POINTER(c_double).in_dll(self._library, 'pr')

    def get_string_buffer(self, length):
        return create_string_buffer(length)

    def space(self, size):
        return self._library.space(size)

    def fold(self, seq, constraints=None, partfunc=False, temperature=None):
        if temperature is not None:
            self.set_temperature(temperature)
        # Structure buffer
        struc_buffer = self.get_string_buffer(len(seq) + 1)
        if constraints is not None:
            self.set_fold_constrained(1)
            struc_buffer.value = constraints
        # Folding
        result = {}
        result['mfe'] = self._library.fold(seq, struc_buffer)
        result['mfe_structure'] = struc_buffer.value
        # Partition function
        if partfunc:
            self._library.init_pf_fold(len(seq))
            if constraints is not None:
                struc_buffer.value = constraints
            else:
                struc_buffer = self.get_string_buffer(len(seq) + 1)
            pffold = self._library.pf_fold(seq, struc_buffer)
            self._library.free_pf_arrays()
            kT = ((temperature + 273.15) * 1.98717) / 1000.
            result['efe_structure'] = struc_buffer.value
            result['efe'] = pffold
            result['mfe_frequency'] = math.exp((pffold - result['mfe']) / kT)
        return result

    def cofold(self, seq1, seq2, constraints=None, partfunc=False, temperature=None):
        if temperature is not None:
            self.set_temperature(temperature)
        # Structure buffer
        struc_buffer = self.get_string_buffer(len(seq1) + len(seq2) + 1)
        if constraints is not None:
            self.set_fold_constrained(1)
            struc_buffer.value = constraints
        # Set cut-point
        cut_point = len(seq1) + 1
        self.set_cut_point(cut_point)
        # Folding
        result = {}
        result['mfe'] = self._library.cofold(seq1+seq2, struc_buffer)
        result['mfe_structure'] = struc_buffer.value[:cut_point-1] + '&' + struc_buffer.value[cut_point-1:]
        # Partition function
        if partfunc:
            self._library.init_co_pf_fold(len(seq1) + len(seq2))
            if constraints is not None:
                struc_buffer.value = constraints
            pffold = self._library.co_pf_fold(seq1+seq2, struc_buffer)
            self._library.free_co_pf_arrays()
            kT = ((temperature + 273.15) * 1.98717) / 1000.
            result['efe_structure'] = struc_buffer.value[:cut_point-1] + '&' + struc_buffer.value[cut_point-1:]
            result['efe'] = pffold.FAB
            result['mfe_frequency'] = math.exp((pffold.FAB - result['mfe']) / kT)
            result['efe_binding'] = pffold.FcAB - pffold.FA - pffold.FB
        return result
