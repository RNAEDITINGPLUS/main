# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Thermodynamics features."""

import math

from . import if_exe_rna
from . import seed

def get_pairing_string(pairing):
    """Returns the pairing with the dots and brackets notation"""
    string = ''
    for i in pairing[::-1]:
        if i == 0:
            string += '.'
        else:
            string += '('
    for i in pairing:
        if i == 0:
            string += '.'
        else:
            string += ')'
    return string

class mmThermo(seed.mmSeed):
    def eval_dg_duplex(self, librna=None, pathrna=None, mirna_start_pairing=None, temperature=None):
        """Computes the *ΔG duplex*, *ΔG binding*, *ΔG seed duplex* and *ΔG seed binding* scores.

           :param librna: Link to the Vienna RNA library.
           :type librna: :class:`LibraryLink`
           :param pathrna: Path to the Vienna RNA executable.
           :type pathrna: str
           :param mirna_start_pairing: Starting position of the seed in the miRNA (from the 5').
           :type mirna_start_pairing: int
           :param temperature: Folding temperature.
           :type temperature: float"""
        # Parameters
        if pathrna is not None:
            if_rna = if_exe_rna.RNAvienna(pathrna)
        elif librna is not None:
            if_rna = librna
        elif hasattr(self, 'libs') and 'rna' in self.libs.libs:
            if_rna = self.libs.get_library_link('rna')
        else:
            if_rna = if_exe_rna.RNAvienna()
        if mirna_start_pairing is None:
            mirna_start_pairing = Defaults.mirna_start_pairing
        if temperature is None:
            temperature = Defaults.temperature
        # Reset
        self.dg_duplex_seeds = []
        self.dg_binding_seeds = []
        self.dg_duplexs = []
        self.dg_duplex_foldings = []
        self.dg_bindings = []
        # Compute
        for its in range(len(self.end_sites)):
            # Target site and seed binding sequences
            target_site_seq = self.target_seq[self.end_sites[its] - self.min_target_length : self.end_sites[its]]
            target_seed_seq = self.target_seq[self.end_sites[its] - (mirna_start_pairing - 1) - self.seed_lengths[its] : self.end_sites[its] - (mirna_start_pairing - 1)]
            mirna_seed_seq = self.mirna_seq[mirna_start_pairing - 1 : mirna_start_pairing + self.seed_lengths[its] - 1]
            # Constraint sequence
            len_no_constraints = self.min_target_length - self.seed_lengths[its] - (mirna_start_pairing - 1)
            constraints_seq = '.' * len_no_constraints + get_pairing_string(self.pairings[its]) + '.' * len_no_constraints
            # Co-folding of seed
            result = if_rna.cofold(target_seed_seq, mirna_seed_seq, partfunc=True, temperature=temperature)
            self.dg_duplex_seeds.append(result['mfe'])
            self.dg_binding_seeds.append(result['efe_binding'])
            # Co-folding of target site
            result = if_rna.cofold(target_site_seq, self.mirna_seq, constraints=constraints_seq, partfunc=True, temperature=temperature)
            self.dg_duplexs.append(result['mfe'])
            self.dg_duplex_foldings.append(result['mfe_structure'])
            self.dg_bindings.append(result['efe_binding'])

    def eval_dg_open(self, librna=None, pathrna=None, upstream_rest=None, downstream_rest=None, dg_binding_area=None, temperature=None):
        """Computes the *ΔG open* score.

           :param librna: Link to the Vienna RNA library.
           :type librna: :class:`LibraryLink`
           :param pathrna: Path to the Vienna RNA executable.
           :type pathrna: str
           :param upstream_rest: Upstream unfolding length.
           :type upstream_rest: int
           :param downstream_rest: Downstream unfolding length.
           :type downstream_rest: int
           :param dg_binding_area: Supplementary sequence length to fold (applied twice: upstream and downstream).
           :type dg_binding_area: int
           :param temperature: Folding temperature.
           :type temperature: float"""
        # Parameters
        if pathrna is not None:
            if_rna = if_exe_rna.RNAvienna(pathrna)
        elif librna is not None:
            if_rna = librna
        elif hasattr(self, 'libs') and 'rna' in self.libs.libs:
            if_rna = self.libs.get_library_link('rna')
        elif hasattr(self, 'exe_path'):
            if_rna = if_exe_rna.RNAvienna(exe_path=self.exe_path)
        else:
            if_rna = if_exe_rna.RNAvienna()
        if upstream_rest is None:
            upstream_rest = Defaults.upstream_rest
        if downstream_rest is None:
            downstream_rest = Defaults.downstream_rest
        if dg_binding_area is None:
            dg_binding_area = Defaults.dg_binding_area
        if temperature is None:
            temperature = Defaults.temperature
        # Reset
        self.dg_opens = []
        # Compute
        for its in range(len(self.end_sites)):
            len_dg_open_seq = dg_binding_area * 2 + upstream_rest + downstream_rest + self.min_target_length
            # Sequence (positions are 1-based ; only when cutting the target sequence, they are changed to 0-based)
            len_polya_upstream = 0
            len_polya_downstream = 0
            start_dg_open_targetseq = 0
            end_dg_open_targetseq = 0
            start_theoretic = self.end_sites[its] - self.min_target_length - upstream_rest - dg_binding_area + 1
            end_theoretic = self.end_sites[its] + downstream_rest + dg_binding_area
            if start_theoretic < 1:
                start_dg_open_targetseq = 1
                len_polya_upstream = abs(start_theoretic) + 1
            else:
                start_dg_open_targetseq = start_theoretic
                len_polya_upstream = 0
            if end_theoretic > self.len_target_seq:
                end_dg_open_targetseq = self.len_target_seq
                len_polya_downstream = end_theoretic - self.len_target_seq
            else:
                end_dg_open_targetseq = end_theoretic
                len_polya_downstream = 0
            seq_for_dg_open = len_polya_upstream * 'A' + self.target_seq[start_dg_open_targetseq - 1 : end_dg_open_targetseq] + len_polya_downstream * 'A'
            # Constraint sequences
            constraints_seq = '.' * dg_binding_area + 'x' * (upstream_rest + self.min_target_length + downstream_rest) + '.' * dg_binding_area
            # Folding
            # dg0
            result_dg0 = if_rna.fold(seq_for_dg_open, partfunc=True, temperature=temperature)
            # dg1
            result_dg1 = if_rna.fold(seq_for_dg_open, constraints=constraints_seq, partfunc=True, temperature=temperature)
            # dg_open
            self.dg_opens.append(result_dg1['efe'] - result_dg0['efe'])

    def eval_dg_total(self):
        """Computes the *ΔG total* score combining *ΔG duplex* and *ΔG open* scores."""
        # Check (and run) dependencies
        if hasattr(self, 'dg_duplexs') is False:
            self.eval_dg_duplex()
        if hasattr(self, 'dg_opens') is False:
            self.eval_dg_open()
        # Reset
        self.dg_totals = []
        # Compute
        for its in range(len(self.end_sites)):
            self.dg_totals.append(self.dg_duplexs[its] + self.dg_opens[its])

    def get_dg_duplex(self, method=None):
        """*ΔG duplex* score with default parameters.

           :param method: Method name used to combine target scores (Example: 'min').
           :type method: str"""
        if method is None:
            method = 'min'
        if hasattr(self, 'dg_duplexs') is False:
            self.eval_dg_duplex()
        if method == 'min':
            self._dg_duplex = min(self.dg_duplexs)
        return self._dg_duplex

    def get_dg_binding(self, method=None):
        """*ΔG binding* score with default parameters.

           :param method: Method name used to combine target scores (Example: 'min').
           :type method: str"""
        if method is None:
            method = 'min'
        if hasattr(self, 'dg_bindings') is False:
            self.eval_dg_duplex()
        if method == 'min':
            self._dg_binding = min(self.dg_bindings)
        return self._dg_binding

    def get_dg_duplex_seed(self, method=None):
        """*ΔG seed duplex* score with default parameters.

           :param method: Method name used to combine target scores (Example: 'min').
           :type method: str"""
        if method is None:
            method = 'min'
        if hasattr(self, 'dg_duplex_seeds') is False:
            self.eval_dg_duplex()
        if method == 'min':
            self._dg_duplex_seed = min(self.dg_duplex_seeds)
        return self._dg_duplex_seed

    def get_dg_binding_seed(self, method=None):
        """*ΔG seed binding* score with default parameters.

           :param method: Method name used to combine target scores (Example: 'min').
           :type method: str"""
        if method is None:
            method = 'min'
        if hasattr(self, 'dg_binding_seeds') is False:
            self.eval_dg_duplex()
        if method == 'min':
            self._dg_binding_seed = min(self.dg_binding_seeds)
        return self._dg_binding_seed

    def get_dg_open(self, method=None):
        """*ΔG open* score with default parameters.

           :param method: Method name used to combine target scores (Example: 'min').
           :type method: str"""
        if method is None:
            method = 'min'
        if hasattr(self, 'dg_opens') is False:
            self.eval_dg_open()
        if method == 'min':
            self._dg_open = min(self.dg_opens)
        elif method == 'exp_sum':
            self._dg_open = 0.
            for dg in self.dg_opens:
                self._dg_open += math.exp(dg * -1.0)
            self._dg_open = math.log(self._dg_open) * -1.0
        return self._dg_open

    def get_dg_total(self, method=None):
        """*ΔG total* score with default parameters.

           :param method: Method name used to combine target scores (Example: 'min').
           :type method: str"""
        if method is None:
            method = 'min'
        if hasattr(self, 'dg_totals') is False:
            self.eval_dg_total()
        if method == 'min':
            self._dg_total = min(self.dg_totals)
        elif method == 'exp_sum':
            self._dg_total = 0.
            for dg in self.dg_totals:
                self._dg_total += math.exp(dg * -1.0)
            self._dg_total = math.log(self._dg_total) * -1.0
        return self._dg_total

    @property
    def dg_duplex(self):
        """*ΔG duplex* score with default parameters."""
        return self.get_dg_duplex()

    @property
    def dg_binding(self):
        """*ΔG binding* score with default parameters."""
        return self.get_dg_binding()

    @property
    def dg_duplex_seed(self):
        """*ΔG seed duplex* score with default parameters."""
        return self.get_dg_duplex_seed()

    @property
    def dg_binding_seed(self):
        """*ΔG seed binding* score with default parameters."""
        return self.get_dg_binding_seed()

    @property
    def dg_open(self):
        """*ΔG open* score with default parameters."""
        return self.get_dg_open()

    @property
    def dg_total(self):
        """*ΔG total* score with default parameters."""
        return self.get_dg_total()

class Defaults(object):
    temperature = 37.
    mirna_start_pairing = 2
    upstream_rest = 10
    downstream_rest = 15
    dg_binding_area = 70
