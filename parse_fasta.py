# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 10:38:41 2016

@author: dippercheng

"""
import os
import re
import numpy as np
import pandas as pd
from Bio import SeqIO

clean = os.listdir('clean')
genus = os.listdir('genus')[1:]
reviewed = os.listdir('reviewed')[1:]


# noinspection PyGlobalUndefined
class Mass:
    def __init__(self, protein_seq):
        self.seq = protein_seq
        self.fragments = [protein_seq]
        self.aa_mass = {
            'A': 71.07884,
            'C': 103.14483,
            'D': 115.08864,
            'E': 129.11552,
            'F': 147.17660,
            'G': 57.05196,
            'H': 137.14120,
            'I': 113.15948,
            'K': 128.17416,
            'L': 113.15948,
            'M': 131.19860,
            'N': 114.10392,
            'P': 97.11672,
            'Q': 128.13080,
            'R': 156.18764,
            'S': 87.07824,
            'T': 101.10512,
            'V': 99.13260,
            'W': 186.21328,
            'Y': 163.17600,
            'X': 110.0,
            'U': 182.08,
            'B': 0.0,
            'O': 255.31,
            'Z': 0.0,
            'J': 0.0
        }

    @property
    def get_mass(self):
        for fragment in self.fragments:
            mass = 18.0153  # - self.aa_mass["M"]
            for aa in fragment:
                mass += self.aa_mass[aa]
        return mass


def parse_description(record):
    key = re.compile(r'\w+\s(.+)\sOS=(\w+)\s(\w+)(.+)GN=(.+)\sPE=[1234]')
    description_pieces = record.description.lower().split()
    locking = key.findall(description_pieces)
    try:
        key_value = locking[0]
        return key_value
    except:
        return ['', '', '', '', 15000, '']


def parse_fasta(species):
    raw_fasta = SeqIO.parse('genus/' + species, "fasta")
    genus_table = pd.DataFrame([parse_description(ref) for ref in raw_fasta],
                               columns=['genus', 'species', 'strain', 'gn', 'mw', 'description'])
    genus_table = genus_table[(genus_table.mw < 12000)]
    genus_table.to_csv("reviewed/" + species[:-6] + ".csv")
    return genus_table


def protein_purge(genus_table_uncleaned):
    genus_groups = [i for i in pd.read_csv('reviewed/' + genus_table_uncleaned).groupby('genus')]
    genus_groups_name = [name[1] for name in genus_groups if genus_table_uncleaned[:-4] == name[0].strip('[]')]
    if len(genus_groups_name) < 4:
        genus_table_cleaned = pd.concat(genus_groups_name)
        genus_table_cleaned['genus'] = genus_table_uncleaned[:-4]
        genus_table_cleaned = genus_table_cleaned[genus_table_cleaned.species != 'sp.'].sort_values(['strain', 'mw'])
        genus_table_cleaned[genus_table_cleaned.pe != 5].to_csv('clean/' + genus_table_uncleaned)
    else:
        print(genus_table_uncleaned + 'failed.')


def get_data(clean_data_name):
    raw_data = pd.read_csv('clean/' + clean_data_name)[['id', 'species', 'gn', 'pe', 'mw']].drop('species' == 'sp.')
    group_species = raw_data.groupby('species')
    step_one = pd.DataFrame([[group[0] for group in group_species],
                             [len(group[1].groupby('gn').median()) for group in group_species]]).T
    step_one.columns = ['species', 'number']
    poly_pecies = step_one[step_one.number > 50].sort_values('number', 0, 0)
    poly_pecies = poly_pecies[poly_pecies.species != 'sp.']
    return poly_pecies



