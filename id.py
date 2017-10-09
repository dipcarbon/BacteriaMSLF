#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 15:22:13 2017

@author: dippercheng
"""

import random
import os
import numpy as np
import pandas as pd
from scipy.stats import mode

'''
groups = random.sample(range(244),244)
group = tuple([groups[i*24:i*24+24] for i in range(6)]+[groups[i*25+24*6-1:i*25+25+24*6-1] for i in range(4)])

group = ([39, 57, 239, 167, 34, 155, 137, 10, 230, 62, 112, 243, 100, 202, 226, 219, 79, 108, 7, 161, 238, 171, 154, 122],\
 [82, 127, 197, 35, 60, 145, 106, 153, 172, 182, 15, 23, 18, 119, 217, 24, 51, 160, 196, 37, 86, 131, 235, 53], \
 [77, 175, 240, 152, 109, 20, 118, 208, 2, 32, 71, 144, 236, 187, 28, 223, 129, 164, 241, 207, 124, 66, 1, 74], \
 [147, 67, 98, 65, 221, 130, 150, 192, 134, 0, 59, 199, 163, 29, 56, 224, 222, 103, 149, 64, 159, 195, 6, 215], \
 [116, 99, 91, 209, 220, 75, 9, 191, 69, 104, 45, 54, 63, 44, 183, 156, 198, 73, 242, 126, 166, 146, 185, 19], \
 [228, 162, 97, 234, 225, 76, 138, 55, 40, 200, 232, 61, 173, 85, 84, 89, 184, 212, 105, 188, 120, 95, 180, 5],\
 [5, 30, 25, 101, 214, 87, 135, 22, 174, 88, 26, 186, 237, 203, 31, 81, 204, 139, 136, 141, 190, 93, 11, 140, 42], \
 [38, 113, 92, 229, 176, 48, 189, 114, 47, 16, 206, 231, 21, 165, 177, 8, 13, 133, 125, 115, 68, 72, 83, 227, 36], \
 [157, 168, 80, 170, 17, 123, 117, 49, 218, 193, 41, 158, 213, 90, 111, 194, 216, 181, 46, 14, 178, 4, 201, 58, 110],\
 [27, 132, 142, 12, 151, 128, 78, 70, 169, 96, 50, 205, 94, 107, 33, 233, 143, 121, 210, 3, 148, 52, 211, 43, 179])
'''
def get_group_sample(validation_id):
    return list(set(range(244))-set(group[validation_id]))

#[224,338,]
# ms_spectra= os.listdir('spectra')[1:]
#lk_spectra = os.listdir('dis_spectra')[1:]
with open('name_spectra.txt') as f:
    ms_spectra = [i.strip() for i in f]
uniprot = os.listdir('uniprot')
the_threshold = 0.12
tolerance = 1000


class BacteriaSpectra:
    def __init__(self, spectra_id):
        self.id = spectra_id
        self.pattern = ms_spectra[spectra_id]
        self.genus = self.pattern.split()[0].lower()
        self.species = self.pattern.split()[1].lower()

    def get_filterd_pattern(self, int_threshold=the_threshold):
        peaks = open('spectra/' + self.pattern, "r")
        standard_mass, standard_int = ([], [])
        for peak in peaks:
            standard_mass.append(float(peak.split()[0]))
            standard_int.append(float(peak.split()[1]))
        peaks.close()
        pattern = pd.DataFrame({'mz': standard_mass, 'int': standard_int})
        if pattern.int.max() > 100:
            pattern['int'] = pattern.int / pattern.int.max()
        elif pattern.int.max() > 1:
            pattern['int'] = pattern.int / 100
        return pattern[pattern.int > int_threshold]

    def get_uniprot_table(self):
        my_table = pd.read_csv('uniprot/' + self.genus + '.csv')[['species', 'gn', 'mw']]
        my_table = my_table[my_table.species == self.species]
        return my_table

    def get_matched_peak(self,int_threshold):
        def match_milk_peak(x):
            y = abs(self.get_filterd_pattern(int_threshold)['mz'] - x)
            y = y[y < x * tolerance * 1e-6]
            if len(y) == 0:
                return np.nan
            return round(y.iloc[0] * 1e6 / x, 2)

        my_table_matched = self.get_uniprot_table()
        my_table_matched['mw'] = my_table_matched['mw'].map(match_milk_peak)
        my_table_matched['genus'] = self.genus
        return my_table_matched.dropna()

    def get_matched_peak_raw(self,int_threshold):
        def match_milk_peak(x):
            raw = self.get_filterd_pattern(int_threshold)['mz']
            y = abs(raw - x)
            y = y[y < x * tolerance * 1e-6]
            if len(y) == 0:
                return np.nan
            return raw[y.index[0]]

        my_table_matched = self.get_uniprot_table()
        my_table_matched['mw'] = my_table_matched['mw'].map(match_milk_peak)
        my_table_matched['genus'] = self.genus
        return my_table_matched.dropna()

def training_gene_names(sample,threshold):
    matched_table = pd.concat([BacteriaSpectra(sample_id).get_matched_peak(threshold) for sample_id in sample],ignore_index=True)
    pre_gn_table = dict(zip(list(set(matched_table.gn)),[0]*len(set(matched_table.gn))))
    for s_iterator in matched_table.groupby(['species','gn']):
        pre_gn_table[s_iterator[0][1]] += np.exp(-np.median(s_iterator[1]['mw'])/1000)

    pre_gn_table = pd.DataFrame(pre_gn_table,index=['value']).T.sort_values('value',ascending=False)
    return (pre_gn_table/max(pre_gn_table.value))


def training_comparison_table(gene_value_table):
    gn_set = list(gene_value_table.index)
    def is_gn_set(x):
        if x in gn_set:
            return x
        return None
    def if_mode(x):
        return round(min(mode(x).mode),1)
    uniprot_table = pd.concat([pd.read_csv('uniprot/' + i)[['species', 'gn', 'mw','genus']] for i in uniprot],ignore_index=True)
    uniprot_table['gn'] = uniprot_table['gn'].map(is_gn_set)
    uniprot_table = uniprot_table.dropna()
    table_iterator =  uniprot_table.groupby(['genus','species'])
    def each_gn(iterator):
        gene_table = dict(zip(['genus','species'],iterator[0]))
        for gn_iterator in iterator[1].groupby('gn'):
            gene_table[gn_iterator[0]] = if_mode(gn_iterator[1]['mw'])
        return pd.DataFrame(gene_table,index=[0])
    comparison_table = pd.concat([each_gn(i) for i in table_iterator],ignore_index=True).round(0)

    return comparison_table.round(0)


class IdentifySpectra(BacteriaSpectra):
    
    def answer(self,model_data,threshold=the_threshold):

        com_table,answer_table,gn_set,gn_value = model_data 
        pattern = self.get_filterd_pattern(threshold)
        
        def match_soy_peak(x):
            y = abs(pattern['mz'] - x)
            y = y[y < x * tolerance * 1e-6]
            if len(y) == 0:
                return 0.0
            return round(np.exp(-y.iloc[0]*1e3/x), 2)
        
        def making_faster(my_table,my_gene_set):
            my_table = my_table.fillna(0)
            b=set({})
            for i in my_gene_set:
                b = b|set(my_table[i].unique())
            return list(b)

        b = making_faster(com_table,gn_set)
        c = [match_soy_peak(x) for x in b]
        my_mw_dict = dict(zip(b,c))
        
        table = com_table.applymap(lambda x:my_mw_dict[x])
        table.to_csv('a_table.csv')
        id_value = table.dot(gn_value)
        id_value.to_csv('a_can.csv')
        score = id_value.max().value
        if score == 0 :
            return ()
        candidate = id_value[id_value.value == score].index
        can = list(set([answer_table.iloc[c]['genus'] for c in candidate]))
        can.sort()
        
        if len(can) == 0:
            return ''
        return can[0]

        #the_answer['score'] = id_value
        #return the_answer.sort_values('score',ascending=False)[:5]
'''
class IdentifySpectra(BacteriaSpectra):
    def answer(self,com_table,answer_table,gn_set,gn_value):
        pattern = self.get_filterd_pattern()
        def match_soy_peak(x):
            
            if x == np.nan:
                return 0.0
            y = abs(pattern['mz'] - x)
            y = y[y < x * tolerance * 1e-6]
            if len(y) == 0:
                return 0.0
            return round(np.exp(-y.iloc[0]*1e3/x), 2)

        table = com_table.applymap(match_soy_peak).dropna(how='all',subset=gn_set)
        id_value = table.dot(gn_value)
        score = id_value.max().value
        if score == 0 :
            return ()
        candidate = id_value[id_value.value == score].index
        can = tuple(set([answer_table.iloc[c]['genus'] for c in candidate]))
        return can
        #the_answer['score'] = id_value
        #return the_answer.sort_values('score',ascending=False)[:5]
'''
def get_accuracy(validation):
    the_accuracy = 0
    for i in range(len(validation)):
        if validation.iloc[i].key in validation.iloc[i].record:
            the_accuracy += 1
    return round(the_accuracy*100/len(validation),2)


#
def gn_training(my_sample,threshold):
    
    my_gene = training_gene_names(my_sample,threshold)[:10]
    my_comparison = training_comparison_table(my_gene)
    my_gene_set = list(my_gene.index)
    my_table = my_comparison[my_gene_set].fillna(0.0)
    my_comparison.to_csv('a_com.csv')
    my_answer_table = my_comparison[['genus','species']]

    model_data = (my_table,my_answer_table,my_gene_set,my_gene)

    return model_data


def validation(my_test,model_data,threshold):
    my_table,my_answer_table,my_gene_set,my_gene = model_data

    my_validation = {}

    for each_sample_id in my_test:
        sample = IdentifySpectra(each_sample_id)
        model_data = (my_table,my_answer_table,my_gene_set,my_gene)
        my_record = sample.answer(model_data,threshold = threshold)
        my_key = sample.genus
        my_validation[each_sample_id] = (my_record,my_key)
    my_validation = pd.DataFrame(my_validation,index=['record','key']).T
    
    my_accuracy = get_accuracy(my_validation)
    print(my_accuracy)
    #my_validation.to_csv(str(my_id)+':'+str(my_accuracy)+'.csv')
    return my_accuracy,my_validation

'''
with open('filtered.txt') as f:
    w = [int(i) for i in f.read().split()]
'''

def int_training(my_sample, my_test, threshold = 0.04):
    #my_sample = group[my_int_id]
    #my_test = group[my_int_idc]

    #my_sample = random.sample(w,10)
    #my_test = my_sample

    model_para = gn_training(my_sample,0.1)
    model_record = validation(my_test,model_para, threshold)

    return model_para,model_record

'''
def th_training(my_id):

    def get_less_group_sample(validation_id):
        return list(set(range(244))-set(group[validation_id])-set(group[(my_id+1)%10]))

    my_sample = group[my_id]
    my_int_sample = get_less_group_sample(my_id)
    my_test = group[(my_id+1)%10]
    my_gene = training_gene_names(my_sample)[:10]
    my_comparison = training_comparison_table(my_gene)
    
    my_gene_set = list(my_gene.index)
    my_table = my_comparison[my_gene_set].fillna(0.0)
    my_answer_table = my_comparison[['genus','species']]
    
    my_validation = {}
    
    for each_sample_id in my_test:
        sample = IdentifySpectra(each_sample_id)
        my_record = ';'.join(sample.answer(my_table,my_answer_table,my_gene_set,my_gene))
        my_key = sample.genus
        my_validation[each_sample_id] = (my_record,my_key)
    my_validation = pd.DataFrame(my_validation,index=['record','key']).T
    
    my_accuracy = get_accuracy(my_validation
    print(my_accuracy)
    #my_validation.to_csv(str(my_id)+':'+str(my_accuracy)+'.csv')
    return (my_gene,my_validation) 


    def get_th_ac(threshold):
        my_validation = {}
        for each_sample_id in th:
            sample = IdentifySpectra(each_sample_id)
            my_record = ';'.join(sample.answer(my_table,my_answer_table,my_gene_set,my_gene,threshold))
            my_key = sample.genus
            my_validation[each_sample_id] = (my_record,my_key)
        my_validation = pd.DataFrame(my_validation,index=['record','key']).T
        my_accuracy = get_accuracy(my_validation)
        return my_accuracy


my = random.sample(w,430)
group = [my[i*43:i*43+43] for i in range(10)]
record = {}
for g in range(4):
    left = [group[gs] for gs in range(10) if gs != g]
    for v in range(9):
        right = random.sample(list(set(my) - set(group[g]) - set(left[v])),50)
        c = int_training(right,left[v])
        #c = int_training([1,2,4], [4,5,6])
        record[str(g) + str(v)] = c
        c[0][-1].to_csv(str(g) + str(v)+'.csv')
        c[1][-1].to_csv(str(g) + str(v)+'r.csv')
'''