# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Probability based on exact distribution feature."""

from . import if_exe_spatt
from . import seed
from . import prob
from . import utils

class mmProbExact(seed.mmSeed):
    def eval_prob_exact(self, libspatt=None, pathspatt=None, markov_order=None, alphabet=None, transitions=None, motif_def=None, motif_upstream_extension=None, motif_downstream_extension=None):
        """Computes the *P.over binomial* score.

           :param libspatt: Link to the Spatt library.
           :type libspatt: :class:`LibraryLink`
           :param pathspatt: Path to the Spatt executable.
           :type pathspatt: str
           :param markov_order: Markov Chain order
           :type markov_order: int
           :param alphabet: List of nucleotides to consider in the sequences (others get filtered).
           :type alphabet: list
           :param transitions: Transition matrix of the Markov Chain model.
           :type transitions: list
           :param motif_def: 'seed' or 'seed_extended' or 'site'.
           :type motif_def: str
           :param motif_upstream_extension: Upstream extension length.
           :type motif_upstream_extension: int
           :param motif_downstream_extension: Downstream extension length.
           :type motif_downstream_extension: int"""
        # Parameters
        if pathspatt is not None:
            if_spatt = if_exe_spatt.Spatt(pathspatt)
        elif libspatt is not None:
            if_spatt = libspatt
        elif hasattr(self, 'libs') and 'spatt' in self.libs.libs:
            if_spatt = self.libs.get_library_link('spatt')
        elif hasattr(self, 'exe_path'):
            if_spatt = if_exe_spatt.Spatt(exe_path=self.exe_path)
        else:
            if_spatt = if_exe_spatt.Spatt()
        if markov_order is None:
            markov_order = Defaults.markov_order
        if alphabet is None:
            try:
                alphabet = Defaults.alphabet
            except AttributeError:
                alphabet = list(set(self.target_seq))
        if transitions is None:
            transitions = prob.get_transitions(self.target_seq, alphabet, markov_order)
        if motif_upstream_extension is None:
            motif_upstream_extension = 0
        if motif_downstream_extension is None:
            motif_downstream_extension = 0
        # Reset
        self.prob_exacts = []
        # Compute
        for its in range(len(self.end_sites)):
            end_site = self.end_sites[its]
            # start_motif and end_motif are sequence coordinates => 1-based
            start_motif, end_motif = seed.get_motif_coordinates(end_site, motif_def, self.pairings[its], motif_upstream_extension, motif_downstream_extension, self.min_target_length)
            motif = self.target_seq[start_motif-1:end_motif]
            self.prob_exacts.append(if_spatt.get_exact_prob(utils.clean_seq(motif, alphabet), self.target_seq.count(motif), self.len_target_seq, alphabet, transitions, markov_order, 'o'))

    def get_prob_exact(self, method=None):
        """*P.over exact* score with default parameters.

           :param method: Method name used to combine target scores (Example: 'min').
           :type method: str"""
        if method is None:
            method = 'min'
        if hasattr(self, 'prob_exacts') is False:
            self.eval_prob_exact()
        if method == 'min':
            self._prob_exact = min(self.prob_exacts)
        elif method == 'prod':
            self._prob_exact = 1.
            for p in self.prob_exacts:
                self._prob_exact *= p
        return self._prob_exact

    @property
    def prob_exact(self):
        """*P.over exact* score with default parameters."""
        return self.get_prob_exact()

class Defaults(object):
    markov_order = 1
