# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Evolutionary features."""

import collections
import copy

import dendropy

from . import if_exe_phast
from . import seed
from . import utils

def get_coord_vec(seq, alphabet, shift=None):
    # By default, coordinates are 1-based
    if shift is None:
        shift = 1
    #coord_vec = array.array('I', [0]*(len(seq)-seq.count('-')))
    coord_vec = []
    i_s = 0
    #i_c = 0
    for s in seq:
        if s in alphabet:
            #coord_vec[i_c] = i_s + 1
            #i_c += 1
            coord_vec.append(i_s + shift)
        i_s += 1
    return coord_vec

def find_all(s, sub, indices=None, offset=None):
    if indices is None:
        indices = []
    if offset is None:
        offset = 0
    i = s.find(sub, offset)
    while i >= 0:
        indices.append(i)
        i = s.find(sub, i + 1)
    return indices

def remove_gap_column(aln):
    clean_aln = copy.copy(aln)
    keep_cols = []
    for i in range(len(aln.values()[0])):
        only_gap = True
        for seq in aln.values():
            if seq[i] != '-':
                only_gap = False
                continue
        if only_gap is False:
            keep_cols.append(i)
    for seq_name in aln.keys():
        clean_aln[seq_name] = ''.join([aln[seq_name][i] for i in keep_cols])
    return clean_aln

class mmEvolution(seed.mmSeed):
    def eval_cons_bls(self, aln_fname=None, aln=None, aln_format=None, aln_alphabet=None, subst_model=None, tree=None, fitting_tree=None, use_em=None, libphast=None, pathphast=None, motif_def=None, motif_upstream_extension=None, motif_downstream_extension=None):
        """Computes the Branch Length Score (*BLS*).

           :param aln_fname: Alignment filename.
           :type aln_fname: str
           :param aln: Alignment it-self.
           :type aln: str
           :param aln_format: Alignment format. Currently supported is FASTA.
           :type aln_format: str
           :param aln_alphabet: List of nucleotides to consider in the aligned sequences (others get filtered).
           :type aln_alphabet: list
           :param subst_model: PhyloFit substitution model (REV...).
           :type subst_model: str
           :param tree: Tree in the Newick format.
           :type tree: str
           :param fitting_tree: Fitting or not the tree on the alignment.
           :type fitting_tree: bool
           :param use_em: Fitting or not the tree with Expectation-Maximization algorithm.
           :type use_em: bool
           :param libphast: Link to the Phast library.
           :type libphast: :class:`LibraryLink`
           :param pathphast: Path to the PHAST executable.
           :type pathphast: str
           :param motif_def: 'seed' or 'seed_extended' or 'site'.
           :type motif_def: str
           :param motif_upstream_extension: Upstream extension length.
           :type motif_upstream_extension: int
           :param motif_downstream_extension: Downstream extension length.
           :type motif_downstream_extension: int"""
        # Parameters
        if aln_fname is None and aln is None:
            raise IOError('An alignment is required')
        if tree is None and fitting_tree is None:
            raise IOError('A tree is required')
        if aln_format is None and aln_fname is not None:
            aln_format = aln_fname.split('.')[-1].upper()
            if aln_format == 'FA':
                aln_format = 'FASTA'
        if aln_format is None:
            raise Exception('Alignment format undetected')
        if aln_alphabet is None:
            aln_alphabet = Defaults.aln_alphabet
        if pathphast is not None:
            if_phast = if_exe_phast.Phast(pathphast)
        elif libphast is not None:
            if_phast = libphast
        elif hasattr(self, 'libs') and 'phast' in self.libs.libs:
            if_phast = self.libs.get_library_link('phast')
        else:
            if_phast = if_exe_phast.Phast()
        if subst_model is None:
            subst_model = Defaults.subst_model
        if fitting_tree is None:
            fitting_tree = True
        if use_em is None:
            use_em = Defaults.use_em
        if motif_upstream_extension is None:
            motif_upstream_extension = 0
        if motif_downstream_extension is None:
            motif_downstream_extension = 0
        fitting_tree_done = False
        # Load alignment
        if aln_fname:
            seqs = utils.load_fasta(aln_fname, upper=True)
        else:
            seqs = utils.load_fasta(aln, as_string=True, upper=True)
        seqs_cleaned = {}
        for seq_name, seq in seqs.items():
            seqs_cleaned[seq_name] = utils.clean_seq(seq, aln_alphabet)
        # Reset
        self.cons_blss = []
        # Compute
        for its in range(len(self.end_sites)):
            end_site = self.end_sites[its]
            # Motif
            # start_motif and end_motif are sequence coordinates => 1-based
            start_motif, end_motif = seed.get_motif_coordinates(end_site, motif_def, self.pairings[its], motif_upstream_extension, motif_downstream_extension, self.min_target_length)
            motif = self.target_seq[start_motif-1:end_motif].replace('U', 'T')
            # Species with seed(s)
            species_with_seed = []
            for seq_name, seq in seqs_cleaned.items():
                if seq.find(motif) != -1:
                    species_with_seed.append(seq_name)
            # BLS
            if len(species_with_seed) > 1:
                # Fitting tree if necessary
                if fitting_tree_done is False:
                    if fitting_tree:
                        if aln_fname:
                            fitted_tree = if_phast.phylofit(subst_model=subst_model, aln_fname=aln_fname, aln_format=aln_format, tree=tree, use_em=use_em)['tree']
                        elif aln:
                            fitted_tree = if_phast.phylofit(subst_model=subst_model, aln=aln, aln_format=aln_format, tree=tree, use_em=use_em)['tree']
                    else:
                        fitted_tree = tree
                    fitting_tree_done = True
                # Compute BLS
                dtree = dendropy.Tree.get_from_string(fitted_tree, schema='newick', preserve_underscores=True)
                dtree.retain_taxa_with_labels(species_with_seed)
                self.cons_blss.append(sum([edge.length for edge in dtree.postorder_edge_iter()][:-1]))
            else:
                self.cons_blss.append(0.)

    def eval_selec_phylop(self, aln_fname=None, aln=None, aln_format=None, aln_alphabet=None, aln_quality=None, mod_fname=None, libphast=None, pathphast=None, method=None, mode=None, motif_def=None, motif_upstream_extension=None, motif_downstream_extension=None):
        """Computes the *PhyloP* score.

           :param aln_fname: Alignment filename.
           :type aln_fname: str
           :param aln: Alignment it-self.
           :type aln: str
           :param aln_format: Alignment format. Currently supported is FASTA.
           :type aln_format: str
           :param aln_alphabet: List of nucleotides to consider in the aligned sequences (others get filtered).
           :type aln_alphabet: list
           :param aln_quality: Check alignment quality (must return True if alignment is fine).
           :type aln_quality: function
           :param mod_fname: Model filename.
           :type mod_fname: str
           :param libphast: Link to the Phast library.
           :type libphast: :class:`LibraryLink`
           :param pathphast: Path to the PHAST executable.
           :type pathphast: str
           :param method: Test name performed by PhyloP (SPH...).
           :type method: str
           :param mode: Testing for conservation (CON), acceleration (ACC) or both (CONACC).
           :type mode: str
           :param motif_def: 'seed' or 'seed_extended' or 'site'.
           :type motif_def: str
           :param motif_upstream_extension: Upstream extension length.
           :type motif_upstream_extension: int
           :param motif_downstream_extension: Downstream extension length.
           :type motif_downstream_extension: int"""
        # Parameters
        if aln_fname is None and aln is None:
            raise IOError('An alignment is required')
        if aln_format is None and aln_fname is not None:
            aln_format = aln_fname.split('.')[-1].upper()
            if aln_format == 'FA':
                aln_format = 'FASTA'
        if aln_format is None:
            raise Exception('Alignment format undetected')
        if aln_alphabet is None:
            aln_alphabet = Defaults.aln_alphabet
        if aln_quality is None:
            aln_quality = Defaults.aln_quality
        if pathphast is not None:
            if_phast = if_exe_phast.Phast(pathphast)
        elif libphast is not None:
            if_phast = libphast
        elif hasattr(self, 'libs') and 'phast' in self.libs.libs:
            if_phast = self.libs.get_library_link('phast')
        elif hasattr(self, 'exe_path'):
            if_phast = if_exe_phast.Phast(exe_path=self.exe_path)
        else:
            if_phast = if_exe_phast.Phast()
        if method is None:
            method = Defaults.method
        if mode is None:
            mode = Defaults.mode
        if motif_upstream_extension is None:
            motif_upstream_extension = 0
        if motif_downstream_extension is None:
            motif_downstream_extension = 0
        # Load alignment and remove gaps
        if aln_fname:
            seqs = utils.load_fasta(aln_fname, upper=True)
        else:
            seqs = utils.load_fasta(aln, as_string=True, upper=True)
        seqs_cleaned = {}
        seqs_coords = {}
        for seq_name, seq in seqs.items():
            seqs_cleaned[seq_name] = utils.clean_seq(seq, aln_alphabet)
            seqs_coords[seq_name] = get_coord_vec(seq, aln_alphabet)
        ref_species = seqs.keys()[0]
        # Reset
        self.selec_phylops = []
        # Compute
        for its in range(len(self.end_sites)):
            end_site = self.end_sites[its]
            # Motif
            # start_motif and end_motif are sequence coordinates => 1-based
            start_motif, end_motif = seed.get_motif_coordinates(end_site, motif_def, self.pairings[its], motif_upstream_extension, motif_downstream_extension, self.min_target_length)
            motif = self.target_seq[start_motif-1:end_motif].replace('U', 'T')
            # Species with seed(s)
            species_with_seed = []
            for seq_name, seq in seqs_cleaned.items():
                if seq.find(motif) != -1:
                    species_with_seed.append(seq_name)
            # phyloP
            pval = 1.0
            if len(species_with_seed) > 1:
                # Extract alignment
                start_motif_in_aln = seqs_coords[ref_species][start_motif-1]
                end_motif_in_aln = seqs_coords[ref_species][end_motif-1]
                partial_seqs = collections.OrderedDict()
                for seq_name, seq in seqs.items():
                    if seq_name in species_with_seed:
                        partial_seqs[seq_name] = seq[start_motif_in_aln - 1:end_motif_in_aln]
                partial_seqs = remove_gap_column(partial_seqs)
                aln = '\n'.join(['> %s\n%s'%(seq_name, seq) for seq_name, seq in partial_seqs.items()])
                if aln_quality(aln):
                    # Compute p-value
                    pval = if_phast.phylop(method=method, mode=mode, mod_fname=mod_fname, aln=aln, aln_format='FASTA')
            self.selec_phylops.append(pval)

    def get_cons_bls(self, method=None):
        """Branch Length Score (*BLS*) score with default parameters.

           :param method: Method name used to combine target scores (Example: 'max').
           :type method: str"""
        if method is None:
            method = 'max'
        if hasattr(self, 'cons_blss') is False:
            self.eval_cons_bls()
        if method == 'max':
            self._cons_bls = max(self.cons_blss)
        elif method == 'avg':
            self._cons_bls = sum(self.cons_blss) / float(len(self.cons_blss))
        return self._cons_bls

    def get_selec_phylop(self, method=None):
        """*PhyloP* score with default parameters.

           :param method: Method name used to combine target scores (Example: 'min').
           :type method: str"""
        if method is None:
            method = 'min'
        if hasattr(self, 'selec_phylops') is False:
            self.eval_selec_phylop()
        if method == 'min':
            self._selec_phylop = min(self.selec_phylops)
        elif method == 'avg':
            self._selec_phylop = sum(self.selec_phylops) / float(len(self.selec_phylops))
        return self._selec_phylop

    @property
    def cons_bls(self):
        """Branch Length Score (*BLS*) score with default parameters."""
        return self.get_cons_bls()

    @property
    def selec_phylop(self):
        """*PhyloP* score with default parameters."""
        return self.get_selec_phylop()

class Defaults(object):
    # PhyloFit
    subst_model = 'REV'
    use_em = True
    # PhyloP
    method = 'SPH'
    mode = 'CONACC'
    #
    aln_alphabet = ['A', 'T', 'C', 'G', 'N']

    @staticmethod
    def aln_quality(aln):
        return True
