# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Interface classe with the `Spatt <http://www.mi.parisdescartes.fr/~nuel/spatt>`_ executable program."""

import os
import re
import subprocess
import tempfile

class Spatt(object):
    """Interface class for the Spatt program."""
    def __init__(self, exe_path=None):
        if exe_path is None:
            self.exe_path = ''
        else:
            self.exe_path = exe_path

    def get_exact_prob(self, motif, nobs, length_seq, alphabet, transitions, markov_order, direction):
        # Cmd
        cmd = [os.path.join(self.exe_path, 'spatt'), '--format', '%a', '--pattern', motif, '--nobs', str(nobs), '--seqlen', str(length_seq), '--alphabet',  ''.join(alphabet), '-m', str(markov_order)]
        # Parameters
        if direction == 'o':
            cmd.append('--over')
        elif direction == 'u':
            cmd.append('--under')
        # Writing markov model
        markovf = tempfile.NamedTemporaryFile(mode='w')
        markovf.write('\n'.join([' '.join(map(str, l)) for l in transitions]))
        markovf.flush()
        cmd.append('--Markov')
        cmd.append(markovf.name)
        # Executing program
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, cwd=tempfile.gettempdir())
        stdout, stderr = p.communicate()
        markovf.close()
        # Parsing results
        decoded = re.search(r'P\(N>=Nobs\)=(?P<prob>\S+)', stdout)
        return float.fromhex(decoded.groupdict()['prob'])
