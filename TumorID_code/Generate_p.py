# -*- coding: utf-8 -*-
"""
@author  : HYH
@E-mail  : huyuhan28@163.com
@className : generate_p
@describe  : TODO  unrelated individual simulation
"""
import pandas as pd
import numpy as np


def gene_alle_file(alle_path):
    '''

    :param alle_path: allele list file generated using frequency files , txt
    :return: two-dimensional list
    '''

    file = open(alle_path)
    loci_ls = []
    allele_data = []
    for line in file.readlines():
        curLine = line.strip().split("\t")
        loci_ls.append(curLine[0])

        curLine.remove(curLine[0])
        floatLine = list(map(float, curLine))
        allele_data.append(floatLine)

    return loci_ls, allele_data


def UP_simulation(alle_path):
    '''

    :param alle_path: allele list file generated using frequency files , txt
    :return:  # dictionary of randome individuals, eg. {'D3S1358': {'Allele 1': 16.0, 'Allele 2': 16.0}, 'vWA': {'Allele 1': 18.0, 'Allele 2': 19.0}, 'D16S539': {'Allele
    '''
    ind_simulation_num = 1
    loci_ls = gene_alle_file(alle_path)[0]
    allele_data = gene_alle_file(alle_path)[1]
    random_ref_dict = {}
    for loci, alle_ls in zip(loci_ls, allele_data):
        geno_sample = sorted(list(np.random.choice(alle_ls, ind_simulation_num * 2)))
        name = ['Allele 1', 'Allele 2']
        loci_dict = dict(zip(name, geno_sample))
        random_ref_dict[loci] = loci_dict
    return random_ref_dict

def RM_txt_simu(alle_path,RM_path):
    '''

    :param alle_path: allele list file generated using frequency files , txt
    :param RM_path: path + genotyping file of randome individuals
    '''

    refDict = UP_simulation(alle_path)

    data_f = pd.DataFrame(refDict)
    data_f = data_f.T
    data_f = data_f.reset_index()

    data_f.rename(columns={'index': 'Marker'}, inplace=True)
    data_f.insert(loc=0,column='Sample Name',value=16)
    data_f.to_csv(RM_path,index=False,sep='\t')

