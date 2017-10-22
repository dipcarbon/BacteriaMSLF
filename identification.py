#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 15:22:13 2017

@author: dippercheng
"""
import itertools
import random
import os
import numpy as np
import pandas as pd
from scipy.stats import mode
import matplotlib.pyplot as plt
from adjustText import adjust_text

with open('spectra_name.txt') as f:
    ms_spectra = [i.strip() for i in f]

uniprot = os.listdir('uniprot')
the_threshold = 0.12
tolerance = 1000


class BacteriaSpectra:
    def __init__(self, spectra_id):
        self.id = spectra_id
        self.pattern = ms_spectra[spectra_id]
        self.taxonomy = self.pattern.split('.')[0]
        self.genus = self.pattern.split()[0].lower()
        self.species = self.pattern.split()[1].lower()

    def plot(self):
        pattern_file = self.get_filterd_pattern(0)
        pat_matched_file = self.get_matched_peak(0.05).groupby('mw_raw').min()
        pat_matched_file = pat_matched_file[pat_matched_file.mw < 500]
        pat_index = pat_matched_file.index
        pattern_file['color'] = pattern_file['mz'].apply(\
            lambda x:[0.6, 0.6, 0.6, 1] if x in pat_index else [0, 0, 0, 1])
        fig, ax = plt.subplots()

        bars = ax.bar(pattern_file['mz'], pattern_file['int'], width=14, facecolor='k')
        
        texts=[]
        for j,rect in enumerate(bars):
            each_mz = pattern_file.iloc[j]['mz']
            if each_mz in pat_index:
                each_text = pat_matched_file.loc[each_mz]['gn']
                left = rect.get_x()+1
                top = rect.get_y()+rect.get_height()+0.01
                texts.append(ax.text(left,top,each_text))
                rect.set_facecolor('#930000')

        adjust_text(texts, add_objects=bars,
                 autoalign='y', only_move={'points':'y', 'text':'y', 'objects':'y'},
                  force_text=0.9, )
        ax.set_ylabel('Intensity')
        ax.set_xlabel('m/z')
        plt.xlim(3000,12000)
        ax.plot([3000,12000],[0,0],color='black')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_title(self.genus.capitalize() +' ' + self.species)
        return ax

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
            raw = self.get_filterd_pattern(int_threshold)
            y = abs(raw['mz'] - x)
            y_delta = y[y < x * tolerance * 1e-6]
            if len(y_delta) == 0:
                return np.nan, np.nan, np.nan
            return round(y_delta.iloc[0] * 1e6 / x, 2), raw.loc[y_delta.index[0]]['mz'], raw.loc[y_delta.index[0]]['int']

        my_table_matched = self.get_uniprot_table()
        ms_interim = my_table_matched['mw'].map(match_milk_peak)
        my_table_matched['mw'] = ms_interim.apply(lambda x:x[0])
        my_table_matched['mw_raw'] = ms_interim.apply(lambda x:x[1])
        my_table_matched['int'] = ms_interim.apply(lambda x:x[2])
        my_table_matched['genus'] = self.genus
        my_table_matched = my_table_matched.dropna()
        my_table_matched.to_csv('matched_result/'+str(self.id)+'.csv')
        return my_table_matched



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
        #table.to_csv('a_table.csv')
        id_value = table.dot(gn_value)
        #id_value.to_csv('a_can.csv')
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
    #my_comparison.to_csv('a_com.csv')
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


def int_training(my_sample, my_test, threshold = 0.04):
    #my_sample = group[my_int_id]
    #my_test = group[my_int_idc]

    #my_sample = random.sample(w,10)
    #my_test = my_sample

    model_para = gn_training(my_sample,0.1)
    model_record = validation(my_test,model_para, threshold)

    return model_para,model_record
#spectra fitlering
'''
with open('spectra_name.txt') as f:
    w = [i for i in f.read().split('\n')]
with open('database_genera.txt') as f:
    genera_list = [i.lower() for i in f.read().split('\n')]

spectra_info = pd.DataFrame(w, columns = ['spectra'])
spectra_info['species'] = spectra_info['spectra'].apply(lambda x: x.split()[1])
spectra_info['genus'] = spectra_info['spectra'].apply(lambda x: x.split()[0].lower())

def if_species_func(x):
    genus = x.split()[0].lower()
    species = x.split()[1]
    if genus not in genera_list:
        return 0
    idid = ms_spectra.index(x)
    table = BacteriaSpectra(idid).get_uniprot_table()
    if len(table) == 0:
        return 1
    else:
        return 2
spectra_info['if_database'] = spectra_info['spectra'].map(if_species_func)
data = spectra_info[spectra_info.if_database == 2]
data ['matched'] = data.index.map(lambda x:len(BacteriaSpectra(x).get_matched_peak(0.05)))
data ['matched'] = data.index.map(lambda x: len(pd.read_csv('matched_result/'+str(x)+'.csv')))
'''
#training demo
'''
sdata = data[data.matched>10]
sdata_list = tuple(sdata.index)

model = gn_training(sdata_list, 0.05)
reslut = validation(sdata_list, model, 0.05)

model = gn_training(random.sample(sdata_list, 30), 0.05)

reslut = validation(abc, model, 0.05)
reslut = validation(sdata_list, model, 0.05)
reslut = validation(random.sample(sdata_list, 30), model, 0.05)
'''


