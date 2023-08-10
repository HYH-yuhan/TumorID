# -*- coding: utf-8 -*-
"""
@author  : HYH
@E-mail  : huyuhan28@163.com
@className : initial_value
@describe  : TODO Calculate the initial value of the parameters with evidence profile
"""
import numpy as np
import pandas as pd
import math

def G_file_process(path):
    '''

    :param path: evidence profile
    :return: dataframe only containing peak height column
    '''
    G_file = pd.read_table('%s' % path)

    # Need to remove gender-related loci
    # corresponding kit :  amle yindel DYS391
    G_file = G_file.drop(index=[5, 6, 10])
    G_file = G_file.drop(
        columns=['Sample Name', 'Marker', 'Allele 1', 'Allele 2', 'Allele 3', 'Allele 4', 'Allele 5', 'Allele 6',
                 'Allele 7', 'Allele 8', 'Allele 9', 'Allele 10', 'Allele 11', 'Allele 12'])
    return G_file

def ini_val(path):
    '''

    :param path: evidence profile
    :return: Expectation , coefficient of variation
    '''

    G_file = G_file_process(path)
    G_file.fillna(0, inplace=True)
    ph_ls = G_file.values.tolist()

    ph_ls = sum(ph_ls, [])

    ph_ls_np = np.array(ph_ls)
    ph_ls_np = ph_ls_np[ph_ls_np != 0]


    E_H = ph_ls_np.mean()
    CV_H = ph_ls_np.std() / E_H

    return E_H,CV_H

def noc_val(path):
    '''

    :param path: evidence profile
    :return: noc
    '''
    G_file = G_file_process(path)
    ph_ls = G_file.values.tolist()

    length = []
    for loci_ph in ph_ls:
        filter_ph_ls = [x for x in loci_ph if not math.isnan(x)]
        length.append(len(filter_ph_ls))

    num = max(length)
    if num < 5 :
        result = 2
    else :
        if (int(num) & 1) == 0:
            result = int(num / 2)
        else:
            result = int((num + 1) / 2)
    return result