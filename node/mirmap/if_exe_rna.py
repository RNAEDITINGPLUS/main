# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Interface classes with the `Vienna RNA <http://www.tbi.univie.ac.at/RNA>`_ executable programs."""

import os
import re
import subprocess
import tempfile

class RNAvienna(object):
    """Interface class for RNA programs from Vienna."""
    def __init__(self, exe_path=None):
        if exe_path is None:
            self.exe_path = ''
        else:
            self.exe_path = exe_path

    def fold(self, seq, constraints=None, partfunc=False, temperature=None):
        return self._fold([seq], 'RNAfold', constraints, partfunc, temperature)

    def cofold(self, seq1, seq2, constraints=None, partfunc=False, temperature=None):
        return self._fold([seq1, seq2], 'RNAcofold', constraints, partfunc, temperature)

    def _fold(self, seqs, prog, constraints=None, partfunc=False, temperature=None):
        # Cmd & regex
        cmd = [os.path.join(self.exe_path, prog), '--noPS']
        regex = r'.+\n(?P<mfe_structure>\S+) \((?P<mfe>.+)\)'
        # Options
        if constraints is not None:
            cmd.append('--constraint')
        else:
            constraints = ''
        if partfunc is True:
            cmd.append('--partfunc')
            if prog == 'RNAfold':
                regex += r'\n *(?P<efe_structure>\S+) \[(?P<efe>.+)\]\n(?P<cfe_structure>\S+) \{ *(?P<cfe>\S+) d=(?P<dist>\S+)\}\n frequency of mfe structure in ensemble (?P<mfe_frequency>\S+); ensemble diversity (?P<ensemble_diversity>\S+)'
            if prog == 'RNAcofold':
                regex += r'\n *(?P<efe_structure>\S+) \[(?P<efe>.+)\]\n frequency of mfe structure in ensemble (?P<mfe_frequency>\S+) , delta G binding= *(?P<efe_binding>\S+)'
        if temperature is not None:
            cmd.append('--temp='+str(temperature))
        # Executing program
        p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, cwd=tempfile.gettempdir())
        stdout, stderr = p.communicate('\n'.join(['&'.join(seqs), constraints]))
        # Parsing results
        decoded = re.match(regex, stdout)
        result = {}
        for k,v in decoded.groupdict().items():
            if 'structure' in k:
                result[k] = v
            else:
                result[k] = float(v)
        return result
