# -*- coding: utf-8 -*-
"""
@author  : HYH
@E-mail  : huyuhan28@163.com
@className : info_read
@describe  : TODO read file
"""
# -*- coding: utf-8 -*-


import pandas as pd

def read_all_fre(path):
    '''
    Read the frequency table
    :param path: path , frequency table
    :return: frequency table dataframe
    '''
    alle_fre_file = pd.read_csv('%s'%path,index_col=0)
    return alle_fre_file


def str_loci_infre(path):
    '''
    Read the loci table
    :param path: path , allele table
    :return:  loci list (str)
    '''
    alle_fre_file = pd.read_csv('%s'%path,index_col=0)
    str_loci_infre = alle_fre_file.columns.tolist()
    return str_loci_infre


def read_sus(path,name):
    '''
    :param path: path ,reference profile
    :param name: reference name 'ref544'
    :return: reference profile dictionary
    '''
    ref_file = pd.read_csv(path)
    ref_file.set_index('SampleName', inplace=True)

    suspect_1 = ref_file.groupby(ref_file.index).get_group(name)

    suspect_1.set_index('Marker', inplace=True)
    susp1 = suspect_1.T.to_dict()
    return susp1

def read_size(sizefile,locus,alle):
    '''

    :param sizefile: size dataframe
    :param locus: locus str
    :param alle:  allele str
    :return: size of specific allele at specific locus
    '''

    df_marker = sizefile.loc[sizefile['Marker'] == locus]
    if alle != 'Q':
        size = float(df_marker.loc[df_marker['Allele'] == float(alle)]['Size'])
    else:
        size = list(df_marker['Marker.Max'])[0]
    return size
