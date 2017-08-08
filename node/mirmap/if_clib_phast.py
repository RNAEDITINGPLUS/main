# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Ctypes interface classes with the `PHAST <http://compgen.bscb.cornell.edu/phast>`_ C library."""

import os

from ctypes import *

from . import if_clib_c

# Classes for structures

# - category_map.h

class CATEGORYMAP(Structure):
    pass

# - lists.h

class LIST(Structure):
    _fields_ = [('array', POINTER(c_void_p)), ('lidx', c_int), ('ridx', c_int), ('CAPACITY', c_int), ('elementsz', c_int), ('step', c_int)]

# - stringplus.h

class STRING(Structure):
    pass

# - gff.h

class GFF_SET(Structure):
    _fields_ = [('features', POINTER(LIST)), ('gff_version', POINTER(STRING)), ('source', POINTER(STRING)), ('source_version', POINTER(STRING)), ('date', POINTER(STRING)), ('groups', POINTER(LIST)), ('group_tag', POINTER(STRING))]

# - list_of_lists.h

class LISTOFLISTS(Structure):
    _fields_ = [('lst', POINTER(LIST)), ('lstName', POINTER(LIST)), ('lstType', POINTER(LIST)), ('class', c_char_p)]

# - msa.h

class MSA(Structure):
    pass

# - tree_model.h

class TM_STRUCT(Structure):
    pass

# - trees.h

class TREE_NODE(Structure):
    pass

# - phylo_fit.h

class PHYLOFIT_STRUCT(Structure):
    _fields_ = [('msa', POINTER(MSA)), ('output_fname_root', c_char_p), ('reverse_group_tag', c_char_p), ('root_seqname', c_char_p), ('subtree_name', c_char_p), ('error_fname', c_char_p), ('see_for_help', c_char_p), ('parsimony_cost_fname', c_char_p), ('msa_fname', c_char_p), ('subst_mod', c_int), ('quiet', c_int), ('nratecats', c_int), ('use_em', c_int), ('window_size', c_int), ('window_shift', c_int), ('use_conditionals', c_int), ('precision', c_int), ('likelihood_only', c_int), ('do_bases', c_int), ('do_expected_nsubst', c_int), ('do_expected_nsubst_tot', c_int), ('do_expected_nsubst_col', c_int), ('random_init', c_int), ('estimate_backgd', c_int), ('estimate_scale_only', c_int), ('do_column_probs', c_int), ('nonoverlapping', c_int), ('gaps_as_bases', c_int), ('no_freqs', c_int), ('no_rates', c_int), ('assume_clock', c_int), ('init_parsimony', c_int), ('parsimony_only', c_int), ('no_branchlens', c_int), ('label_categories', c_int), ('symfreq', c_int), ('init_backgd_from_data', c_int), ('use_selection', c_int), ('max_em_its', c_int), ('nsites_threshold', c_uint), ('tree', POINTER(TREE_NODE)), ('cm', POINTER(CATEGORYMAP)), ('nooptstr', POINTER(STRING)), ('cats_to_do_str', POINTER(LIST)), ('window_coords', POINTER(LIST)), ('ignore_branches', POINTER(LIST)), ('alt_mod_str', POINTER(LIST)), ('bound_arg', POINTER(LIST)), ('rate_consts', POINTER(LIST)), ('label_str', POINTER(LIST)), ('label_type', POINTER(LIST)), ('alpha', c_double), ('selection', c_double), ('gff', POINTER(GFF_SET)), ('input_mod', POINTER(TM_STRUCT)), ('logf', POINTER(if_clib_c.FILE)), ('results', POINTER(LISTOFLISTS))]

# - phylo_p.h

class PHYLOP_STRUCT(Structure):
    _fields_ = [('msa', POINTER(MSA)), ('prior_only', c_int), ('post_only', c_int), ('quantiles_only', c_int), ('output_wig', c_int), ('output_gff', c_int), ('nsites', c_int), ('fit_model', c_int), ('base_by_base', c_int), ('default_epsilon', c_int), ('refidx', c_int), ('refidx_feat', c_int), ('ci', c_double), ('epsilon', c_double), ('subtree_name', c_char_p), ('chrom', c_char_p), ('branch_name', POINTER(LIST)), ('feats', POINTER(GFF_SET)), ('method', c_int), ('mode', c_int), ('outfile', POINTER(if_clib_c.FILE)), ('logf', POINTER(if_clib_c.FILE)), ('mod', POINTER(TM_STRUCT)), ('cats_to_do', POINTER(LIST)), ('cm', POINTER(CATEGORYMAP)), ('help', c_char_p), ('mod_fname', c_char_p), ('msa_fname', c_char_p), ('results', POINTER(LISTOFLISTS)), ('no_prune', c_int)]

class Phast(object):
    """Interface class for the Phast library"""
    def __init__(self, library_path=None, library_name=None):
        if library_name is None:
            lib = 'libphast.so'
        else:
            lib = library_name
        if library_path is not None:
            lib = os.path.join(library_path, lib)
        self._library = cdll.LoadLibrary(lib)
        # Lib C interface
        self._libc = if_clib_c.LibC()
        # Enums
        self._enum_method_types = {'SPH':c_int(0), 'LRT':c_int(1), 'SCORE':c_int(2), 'GERP':c_int(3)}
        self._enum_mode_types = {'CON':c_int(0), 'ACC':c_int(1), 'NNEUT':c_int(2), 'CONACC':c_int(3)}
        self._enum_msa_format_types = {'PHYLIP':c_int(0), 'MPM':c_int(1), 'FASTA':c_int(2), 'SS':c_int(3), 'LAV':c_int(4), 'MAF':c_int(5), 'UNKNOWN_FORMAT':c_int(6)}
        self._enum_list_element_type = {'INT_LIST':c_int(0), 'DBL_LIST':c_int(1), 'CHAR_LIST':c_int(2), 'LIST_LIST':c_int(3)}
        # Functions arguments and result types
        # - gff.h
        self._library.gff_free_set.argtypes = [POINTER(GFF_SET)]
        self._library.gff_free_set.restype = None
        self._library.gff_read_set.argtypes = [POINTER(if_clib_c.FILE)]
        self._library.gff_read_set.restype = POINTER(GFF_SET)
        # - lists.h
        self._library.lst_free.argtypes = [POINTER(LIST)]
        self._library.lst_free.restype = None
        self._library.lst_free_strings.argtypes = [POINTER(LIST)]
        self._library.lst_free_strings.restype = None
        self._library.lst_size_non_inline.argtypes = [POINTER(LIST)]
        self._library.lst_size_non_inline.restype = c_int
        self._library.lst_get_ptr_non_inline.argtypes = [POINTER(LIST), c_int]
        self._library.lst_get_ptr_non_inline.restype = c_void_p
        self._library.lst_get_int_non_inline.argtypes = [POINTER(LIST), c_int]
        self._library.lst_get_int_non_inline.restype = c_int
        self._library.lst_get_dbl_non_inline.argtypes = [POINTER(LIST), c_int]
        self._library.lst_get_dbl_non_inline.restype = c_double
        # - list_of_lists.h
        self._library.lol_free.argtypes = [POINTER(LISTOFLISTS)]
        self._library.lol_free.restype = None
        self._library.lol_new.argtypes = [c_int]
        self._library.lol_new.restype = POINTER(LISTOFLISTS)
        self._library.lol_find_list.argtypes = [POINTER(LISTOFLISTS), c_char_p, c_int]
        self._library.lol_find_list.restype = POINTER(LIST)
        self._library.lol_find_lol.argtypes = [POINTER(LISTOFLISTS), c_char_p]
        self._library.lol_find_lol.restype = POINTER(LISTOFLISTS)
        # - maf.h
        self._library.maf_read_cats.argtypes = [POINTER(if_clib_c.FILE), POINTER(if_clib_c.FILE), c_int, c_char_p, POINTER(GFF_SET), POINTER(CATEGORYMAP), c_int, c_int, c_char_p, c_int, c_int, POINTER(LIST)]
        self._library.maf_read_cats.restype = POINTER(MSA)
        # - memory_handler.h
        self._library.phast_free_all.argtypes = None
        self._library.phast_free_all.restype = None
        self._library.phast_new_mem_handler.argtypes = None
        self._library.phast_new_mem_handler.restype = None
        self._library.phast_num_mem_handlers.argtypes = None
        self._library.phast_num_mem_handlers.restype = c_int
        # - misc.h
        self._library.fopen_fname.argtypes = [c_char_p, c_char_p]
        self._library.fopen_fname.restype = POINTER(if_clib_c.FILE)
        self._library.get_arg_list.argtypes = [c_char_p]
        self._library.get_arg_list.restype = POINTER(LIST)
        self._library.set_seed.argtypes = [c_int]
        self._library.set_seed.restype = None
        # - msa.h
        self._library.msa_free.argtypes = [POINTER(MSA)]
        self._library.msa_free.restype = None
        self._library.msa_new_from_file_define_format.argtypes = [POINTER(if_clib_c.FILE), c_int, c_char_p]
        self._library.msa_new_from_file_define_format.restype = POINTER(MSA)
        # - phylo_fit.h
        self._library.phyloFit_struct_new.argtypes = [c_int]
        self._library.phyloFit_struct_new.restype = POINTER(PHYLOFIT_STRUCT)
        self._library.run_phyloFit.argtypes = [POINTER(PHYLOFIT_STRUCT)]
        self._library.run_phyloFit.restype = c_int
        # - phylo_p.h
        self._library.phyloP_struct_new.argtypes = [c_int]
        self._library.phyloP_struct_new.restype = POINTER(PHYLOP_STRUCT)
        self._library.phyloP.argtypes = [POINTER(PHYLOP_STRUCT)]
        self._library.phyloP.restype = None
        # - subst_mods.h
        self._library.tm_get_subst_mod_type.argtypes = [c_char_p]
        self._library.tm_get_subst_mod_type.restype = c_int
        # - tree_model.h
        self._library.tm_free.argtypes = [POINTER(TM_STRUCT)]
        self._library.tm_free.restype = None
        self._library.tm_new_from_file.argtypes = [POINTER(if_clib_c.FILE), c_int]
        self._library.tm_new_from_file.restype = POINTER(TM_STRUCT)
        # - trees.h
        self._library.tr_new_from_string.argtypes = [c_char_p]
        self._library.tr_new_from_string.restype = POINTER(TREE_NODE)

    def get_phylofit_struct(self, rphast):
        return self._library.phyloFit_struct_new(rphast)

    def get_phylop_struct(self, rphast):
        return self._library.phyloP_struct_new(rphast)

    def clean_phylofit(self, f):
        raise NotImplementedError()

    def clean_phylop(self, p):
        if p.contents.subtree_name is not None:
            self._libc.free(p.contents.subtree_name)
        if p.contents.chrom is not None:
            self._libc.free(p.contents.chrom)
        if p.contents.branch_name:
            self._library.lst_free_strings(p.contents.branch_name)
            self._library.lst_free(p.contents.branch_name)
        if p.contents.msa:
            self._library.msa_free(p.contents.msa)
        if p.contents.feats:
            self._library.gff_free_set(p.contents.feats)
        if p.contents.results:
            self._library.lol_free(p.contents.results)
        if p.contents.mod:
            self._library.tm_free(p.contents.mod)
        self._libc.free(p)

    def phylofit(self, subst_model, aln_fname=None, aln=None, aln_format=None, tree=None, use_em=None, use_mem_handler=None):
        if use_mem_handler is None and self._library.phast_num_mem_handlers() != -1:
            use_mem_handler = True
        if use_mem_handler:
            self._library.phast_new_mem_handler()
        f = self.get_phylofit_struct(0)
        self._library.set_seed(-1)
        f.contents.results = self._library.lol_new(2)
        # Tree
        if tree is not None:
            f.contents.tree = self._library.tr_new_from_string(tree)
        # Substitution model
        f.contents.subst_mod = self._library.tm_get_subst_mod_type(subst_model)
        # EM
        if use_em is True:
            f.contents.use_em = 1
        # MSA
        if aln_fname is not None and aln is None:
            f.contents.msa_fname = aln_fname
            aln_file = if_clib_c.RegularFile(aln_fname, 'r')
        elif aln_fname is None and aln is not None:
            aln_file = if_clib_c.InMemoryFile(aln)
        if aln_format == 'MAF':
            raise NotImplementedError()
        else:
            f.contents.msa = self._library.msa_new_from_file_define_format(aln_file.file, self._enum_msa_format_types[aln_format], 'ACGT')
        aln_file.close()
        # Output
        f.contents.output_fname_root = None
        f.contents.quiet = 1
        # Call phyloFit
        self._library.run_phyloFit(f)
        # Output
        result = {}
        trees = self._library.lol_find_list(f.contents.results, 'tree', self._enum_list_element_type['CHAR_LIST'])
        result['tree'] = cast(self._library.lst_get_ptr_non_inline(trees, 0), c_char_p).value
        # Cleaning
        if use_mem_handler:
            self._library.phast_free_all()
        else:
            self.clean_phylofit(f)
        return result

    def phylop(self, method, mode, mod_fname, aln_fname=None, aln=None, aln_format=None, gff_fname=None, branch=None, prune=None, use_mem_handler=None):
        if use_mem_handler is None and self._library.phast_num_mem_handlers() != -1:
            use_mem_handler = True
        if use_mem_handler:
            self._library.phast_new_mem_handler()
        p = self.get_phylop_struct(0)
        p.contents.results = self._library.lol_new(20)
        p.contents.method = self._enum_method_types[method]
        p.contents.mode = self._enum_mode_types[mode]
        p.contents.mod_fname = mod_fname
        p.contents.msa_fname = aln_fname
        # GFF
        if gff_fname:
            gff_f = self._library.fopen_fname(gff_fname, 'r')
            p.contents.feats = self._library.gff_read_set(gff_f)
            self._libc.fclose(gff_f)
        # Model
        model_f = self._library.fopen_fname(p.contents.mod_fname, 'r')
        p.contents.mod = self._library.tm_new_from_file(model_f, 1)
        self._libc.fclose(model_f)
        # MSA
        if aln_fname is not None and aln is None:
            aln_file = if_clib_c.RegularFile(aln_fname, 'r')
        elif aln_fname is None and aln is not None:
            aln_file = if_clib_c.InMemoryFile(aln)
        if aln_format == 'MAF':
            raise NotImplementedError()
            #p.contents.msa = self._library.maf_read_cats(msa_f, None, 1, None, None, None, -1, 1, None, 0, 0, None)
        else:
            p.contents.msa = self._library.msa_new_from_file_define_format(aln_file.file, self._enum_msa_format_types[aln_format], 'ACGT')
        aln_file.close()
        # Branch
        if branch is not None:
            p.branch_name = self._library.get_arg_list(branch)
        # Output
        p.contents.outfile = None
        # Parameters
        if prune is True:
            p.contents.no_prune = 1
        elif prune is False:
            p.contents.no_prune = 0
        # Call phyloP
        self._library.phyloP(p)
        # Output
        pvals = self._library.lol_find_list(p.contents.results, 'pval.cons', self._enum_list_element_type['DBL_LIST'])
        pval = self._library.lst_get_dbl_non_inline(pvals, 0)
        # Cleaning
        if use_mem_handler:
            self._library.phast_free_all()
        else:
            self.clean_phylop(p)
        return pval
