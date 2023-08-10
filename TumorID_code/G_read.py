# -*- coding: utf-8 -*-
"""
@author  : HYH
@E-mail  : huyuhan28@163.com
@className : G_read
@describe  : TODO Read sample profile（.txt） and convert to dictionary,{'D7S820': {'Allele 1': '8', 'Allele 2': '9', 'Allele 3': 11.0, 'Allele 4': *, 'Allele 5': *, ...}...}
"""


import pandas as pd

def g_read(path):
    '''

    :param path: sample profile path
    :return: dictionary :  {'D7S820': {'Allele 1': '8', 'Allele 2': '9', 'Allele 3': 11.0, 'Allele 4': *, 'Allele 5': *, ...}...}
    '''

    G_file = pd.read_table('%s'%path,encoding="utf-8")
    G_file.set_index('Sample Name', inplace=True)
    G_file.set_index('Marker', inplace=True)

    G_file.fillna('*', inplace=True)

    G = G_file.T.to_dict()

    return G

