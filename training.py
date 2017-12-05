#uniprot_table = pd.concat([pd.read_csv('uniprot/' + i)[['species', 'gn', 'mw','genus']] for i in uniprot],ignore_index=True)
import itertools
import random
import time
import os
import numpy as np
import pandas as pd
from scipy.stats import mode
import matplotlib.pyplot as plt
from adjustText import adjust_text

def training_gene_names(sample, threshold):
    matched_table = pd.concat([BacteriaSpectra(sample_id).get_matched_peak(threshold) for sample_id in sample],ignore_index=True)
    pre_gn_table = dict(zip(list(set(matched_table.gn)),[0]*len(set(matched_table.gn))))
    for s_iterator in matched_table.groupby(['species','gn']):
        pre_gn_table[s_iterator[0][1]] += np.exp(-np.median(s_iterator[1]['mw'])/1000)

    pre_gn_table = pd.DataFrame(pre_gn_table,index=['value']).T.sort_values('value',ascending=False)
    return (pre_gn_table / max(pre_gn_table.value))['value']


def genetic_algorithms(each_group, gene_int):
    gene_out  = gene_int
    return gene_out


def training_comparison_table(gene_value_table):
    gn_set = list(gene_value_table.index)
    def is_gn_set(x):
        if x in gn_set:
            return x
        return None
    def if_mode(x):
        return round(min(mode(x).mode),1)
    uniprot_table = pd.read_csv('uniprot_table.csv', index_col = 0)
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


def validation(my_test, model_data, threshold):
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

def gene_to_model(my_gene):
    my_gene = my_gene[:10]
    my_comparison = training_comparison_table(my_gene)
    my_gene_set = list(my_gene.index)
    my_table = my_comparison[my_gene_set].fillna(0.0)
    my_answer_table = my_comparison[['genus','species']]
    model_data = (my_table,my_answer_table,my_gene_set,my_gene)
    return  model_data

def fold_10(lists):
    train = set(lists)
    size = len(lists)//10
    groups = []
    for i in range(9):
        test = set(random.sample(train, size))
        groups.append(test)
        train = train - test
    groups.append(train)
    return tuple(groups)

def double_cross_vali():
    out_test = random.sample(data_index, 40)
    train = set(data_index) - set(test)

def loop():
    train_group = set(data_index) - set(each_group)
    test_group = set(each_group)
    gene_int = training_gene_names(train_group, 0.1)
    gene_out = genetic_algorithms(train_group, gene_int)
    model_data = gene_to_model(gene_out)
    rate, answer_list = validation(test_group, model_data, 0.1)


def cross_vali():
    group = fold_10(data_index)

    def loop():
        train_group = set(data_index) - set(each_group)
        test_group = set(each_group)
        gene_int = training_gene_names(train_group, 0.1)
        gene_out = genetic_algorithms(train_group, gene_int)
        model_data = gene_to_model(gene_out)
        rate, answer_list = validation(test_group, model_data, 0.1)

"""
t1 = time.time()
threshold=the_threshold
com_table,answer_table,gn_set,gn_value = model_data 
pattern = sample.get_filterd_pattern(threshold)
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
#return (table.dot(gn_value))
id_value_raw = (table.dot(gn_value))
id_value = id_value_raw.sort_values()
#id_value.to_csv('a_can.csv')
id_panel = id_value[-4:]
score_mode =id_panel.mode()
if not score_mode.empty :
    score = score_mode.iloc[-1]
else:
    score = id_panel.max()
#return id_value
#if score == 0:
#   return ()
candidates = id_panel[id_panel == score].index[0]
answer_candidate = answer_table.iloc[candidates]['genus']
t2 = time.time()
print(answer_candidate, score, round(t2-t1, 2), 's')
calculated_table = table
my_answer_table = answer_table.copy()
my_answer_table['score'] = id_value_raw

data = pd.read_csv('newdata.csv',index_col=0)[['species', 'genus']]
data_index = set(data.index)
"""




