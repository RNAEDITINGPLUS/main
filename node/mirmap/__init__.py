# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

""":class:`mm` and :class:`mmPP` base classes of :mod:`miRmap` that inherit their methods from all the modules. Each module define the methods for one category."""

from . import model
from . import prob_binomial
from . import report
from . import targetscan

from . import evolution
from . import prob_exact
from . import thermo

class mm(evolution.mmEvolution, model.mmModel, prob_binomial.mmProbBinomial, prob_exact.mmProbExact, report.mmReport, thermo.mmThermo, targetscan.mmTargetScan):
    """miRNA and mRNA containing class."""
    pass

class mmPP(model.mmModel, prob_binomial.mmProbBinomial, report.mmReport, targetscan.mmTargetScan):
    """miRNA and mRNA containing class with pure Python methods only."""
    pass
