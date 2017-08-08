# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Reporting methods."""

from . import seed

def get_mirna_pairing_string(pairing):
    """Returns the pairing with the dots and pipes notation"""
    string = ''
    for i in pairing[::-1]:
        if i == 0:
            string += '.'
        else:
            string += '|'
    return string

class mmReport(seed.mmSeed):
    def report(self):
        """Returns a formatted report of already computed features for all target site(s)."""
        report_lines = []
        for its in range(len(self.end_sites)):
            end_site = self.end_sites[its]
            start = max(0, end_site-self.len_mirna_seq-10)
            end = end_site + 10
            report_lines.append(str(start+1) + ' '*(end_site-start-len(str(start+1))-1) + str(end_site))
            report_lines.append('|' + ' '*(end_site-start-2) + '|')
            report_lines.append(self.target_seq[start:end])
            seed_pairing_string = get_mirna_pairing_string(self.pairings[its])
            report_lines.append(' ' * (end_site - len(seed_pairing_string) - start) + seed_pairing_string)
            report_lines.append(' ' * (end_site - self.len_mirna_seq - start) + self.mirna_seq[::-1])
	    '''
            try:
                report_lines.append('  %-30s% .5f'%('ΔG duplex (kcal/mol)', self.dg_duplexs[its]))
            except AttributeError:
                pass
            try:
                report_lines.append('  %-30s% .5f'%('ΔG binding (kcal/mol)', self.dg_bindings[its]))
            except AttributeError:
                pass
            try:
                report_lines.append('  %-30s% .5f'%('ΔG open (kcal/mol)', self.dg_opens[its]))
            except AttributeError:
                pass
            try:
                report_lines.append('  %-30s% .5f'%('ΔG total (kcal/mol)', self.dg_totals[its]))
            except AttributeError:
                pass
            try:
                report_lines.append('  %-30s% .5f'%('AU content', self.tgs_aus[its]))
            except AttributeError:
                pass
            try:
                report_lines.append('  %-30s% .5f'%('UTR position', self.tgs_positions[its]))
            except AttributeError:
                pass
            try:
                report_lines.append('  %-30s% .5f'%('3\' pairing', self.tgs_pairing3ps[its]))
            except AttributeError:
                pass
            try:
                report_lines.append('  %-30s% .5f'%('TargetScan score', self.tgs_scores[its]))
            except AttributeError:
                pass
            try:
                report_lines.append('  %-30s% .5f'%('Probability (Exact)', self.prob_exacts[its]))
            except AttributeError:
                pass
            try:
                report_lines.append('  %-30s% .5f'%('Probability (Binomial)', self.prob_binomials[its]))
            except AttributeError:
                pass
            try:
                report_lines.append('  %-30s% .5f'%('Conservation (BLS)', self.cons_blss[its]))
            except AttributeError:
                pass
            try:
                report_lines.append('  %-30s% .5f'%('Conservation (PhyloP)', self.selec_phylops[its]))
            except AttributeError:
                pass
            '''
        return '\t'.join(report_lines)
