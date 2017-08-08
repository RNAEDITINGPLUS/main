# -*- coding: utf-8 -*-

#
# Copyright (C) 2011-2012 Charles E. Vejnar
#
# This is free software, licensed under the GNU General Public License v3.
# See /LICENSE for more information.
#

"""Basic functions for probability features: Transitions matrix, motif probability..."""

from . import utils

def permutations(items, n):
    if n == 0:
        yield ''
    else:
        for i in range(len(items)):
            for base in permutations(items, n - 1):
                yield str(items[i]) + str(base)

def get_transitions(seq, alphabet, markov_order):
    """Computes transitions matrix"""
    transitions = []
    motifs = list(permutations(alphabet, markov_order+1))
    counts = dict(zip(motifs, [0]*len(motifs)))
    for i in range(len(seq)-markov_order):
        s = seq[i:i+markov_order+1]
        if all([i in alphabet for i in s]):
            counts[s] += 1
    for motif in motifs:
        transitions.append(counts[motif])
    transitions = list(utils.grouper(len(alphabet), transitions))
    sums = map(sum, transitions)
    for i in range(len(transitions)):
        if sums[i] != 0:
            transitions[i] = [transitions[i][j] / float(sums[i]) for j in range(len(transitions[i]))]
        else:
            transitions[i] = [0.] * len(transitions[i])
    return transitions

def prob_motif(motif, alphabet, markov_order, transitions):
    """Computes the probability of a motif based on a transitions matrix"""
    transitions = utils.flatten(transitions)
    motifs = list(permutations(alphabet, markov_order+1))
    motifs_index = dict([(motifs[i], i) for i in range(len(motifs))])
    prob = 1.
    for i in range(len(motif) - markov_order):
        prob *= transitions[motifs_index[motif[i:i+markov_order+1]]]
    return prob
