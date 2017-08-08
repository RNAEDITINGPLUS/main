# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Interface classes with the `PHAST <http://compgen.bscb.cornell.edu/phast>`_ executable programs."""

import os
import re
import subprocess
import tempfile

class Phast(object):
    """Interface class for RNA programs from Vienna."""
    def __init__(self, exe_path=None):
        if exe_path is None:
            self.exe_path = ''
        else:
            self.exe_path = exe_path

    def phylofit(self, subst_model=None, aln_fname=None, aln=None, aln_format=None, tree=None, use_em=None):
        # Cmd
        cmd = [os.path.join(self.exe_path, 'phyloFit'), '--precision', 'HIGH', '--out-root', '-', '--msa-format', aln_format]
        # Options
        if subst_model is not None:
            cmd.append('--subst-mod')
            cmd.append(subst_model)
        if use_em is True:
            cmd.append('--EM')
        # Tree
        tree_file = None
        if tree is not None:
            tree_file = tempfile.NamedTemporaryFile(mode='w')
            tree_file.write(tree)
            tree_file.flush()
            cmd.append('--tree')
            cmd.append(tree_file.name)
        # MSA
        aln_file = None
        if aln_fname is not None and aln is None:
            cmd.append(aln_fname)
        elif aln_fname is None and aln is not None:
            aln_file = tempfile.NamedTemporaryFile(mode='w')
            aln_file.write(aln)
            aln_file.flush()
            cmd.append(aln_file.name)
        # Executing program
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=tempfile.gettempdir())
        stdout, stderr = p.communicate()
        if tree_file is not None:
            tree_file.close()
        if aln_file is not None:
            aln_file.close()
        # Parsing results
        decoded = re.match(r'ALPHABET: (?P<alphabet>[^\n]+)\nORDER: (?P<order>\S+)\nSUBST_MOD: (?P<subst_mod>\S+)\nTRAINING_LNL: (?P<training_lnl>\S+)\nBACKGROUND: (?P<background>[^\n]+)\nRATE_MAT:\n(?P<rate_mat>.+)\nTREE: (?P<tree>.+;)', stdout, re.DOTALL)
        result =  decoded.groupdict()
        result['mod_raw'] = stdout
        result['training_lnl'] = float(result['training_lnl'])
        return result

    def phylop(self, method, mode, mod_fname, aln_fname=None, aln=None, aln_format=None, gff_fname=None, branch=None, prune=None):
        # Cmd
        cmd = [os.path.join(self.exe_path, 'phyloP'), '--method', method, '--mode', mode, '--msa-format', aln_format]
        # Options
        if gff_fname is not None:
            cmd.append('--features')
            cmd.append(gff_fname)
        if branch is not None:
            cmd.append('--branch')
            cmd.append(branch)
        if prune is False:
            cmd.append('--no-prune')
        # Model
        cmd.append(mod_fname)
        aln_file = None
        # MSA
        if aln_fname is not None and aln is None:
            cmd.append(aln_fname)
        elif aln_fname is None and aln is not None:
            aln_file = tempfile.NamedTemporaryFile(mode='w')
            aln_file.write(aln)
            aln_file.flush()
            cmd.append(aln_file.name)
        # Executing programs
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=tempfile.gettempdir())
        stdout, stderr = p.communicate()
        if aln_file is not None:
            aln_file.close()
        # Parsing results
        decoded = re.search(r'p-value of conservation: (?P<prob>\S+)', stdout)
        return float(decoded.groupdict()['prob'])
