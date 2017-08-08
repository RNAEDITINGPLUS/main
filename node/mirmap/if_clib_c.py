# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Module collecting interface classes with C library."""

from ctypes import *
from ctypes.util import find_library

# Classes for structures
class FILE(Structure):
    pass

class LibC(object):
    """Interface class for the Lib C"""
    def __init__(self, location=None):
        if location is None:
            self._library = cdll.LoadLibrary(find_library('c'))
        else:
            self._library = cdll.LoadLibrary(location)
        # Functions arguments and result types
        self._library.free.argtypes = [c_void_p]
        self._library.free.restype = None
        self._library.fopen.argtypes = [c_char_p, c_char_p]
        self._library.fopen.restype = POINTER(FILE)
        self._library.fclose.argtypes = [POINTER(FILE)]
        self._library.fclose.restype = c_int
        self._library.open_memstream.argtypes = [POINTER(c_char_p), POINTER(c_size_t)]
        self._library.open_memstream.restype = POINTER(FILE)

    def free(self, pointer):
        self._library.free(pointer)

    def fopen(self, path, mode):
        return self._library.fopen(path, mode)

    def fclose(self, fp):
        return self._library.fclose(fp)

    def fprintf(self, stream, formatstring, *strings):
        return self._library.fprintf(stream, formatstring, *strings)

    def open_memstream(self, ptr, sizeloc):
        return self._library.open_memstream(ptr, sizeloc)

class RegularFile(object):
    """Regular file"""
    def __init__(self, fname, mode, libc=None):
        if libc is None:
            self._libc = LibC()
        else:
            self._libc = libc
        self._file = self._libc.fopen(fname, mode)

    def close(self):
        self._libc.fclose(self._file)

    @property
    def file(self):
        return self._file

class InMemoryFile(object):
    """In-memory file"""
    def __init__(self, stream, libc=None):
        if libc is None:
            self._libc = LibC()
        else:
            self._libc = libc
        self._buf = c_char_p()
        self._buf_len = c_size_t()
        self._file = self._libc.open_memstream(byref(self._buf), byref(self._buf_len))
        self._libc.fprintf(self._file, '%s', stream)

    def close(self):
        self._libc.fclose(self._file)
        self._libc.free(self._buf)

    @property
    def file(self):
        return self._file
