# -*- coding: utf-8 -*-
"""
@author  : HYH
@E-mail  : huyuhan28@163.com
@className : G_read_allePH
@describe  : TODO adjust the dictionary format:  {'D7S820': {8.0: 499.0, 9.0: 105.0, 11.0: 246.0}  ...}
"""
# -*- coding: utf-8 -*-


import G_read
import logging


def g_read_allePH(G):
    '''

    :param G:   the  dictionary before adjustment , Generated by G_read
    :return:    the  adjusted  dictionary , {'D7S820': {8.0: 499.0, 9.0: 105.0, 11.0: 246.0}  ...}
    '''

    strloci_g = list(G.keys())

    alle_peakHeight_list =[]
    for strindex in G.values():   #  each locus

        #  create allele list
        alle_list = []
        height_list = []

        for i in range(1, 12):

            if strindex['Allele %s' % i] != '*':
                alle_list.append(strindex['Allele %s' % i])
            if strindex['Height %s' % i] != '*':
                height_list.append(strindex['Height %s' % i])

        #  convert number to float
        height_list = list(map(float, height_list))
        try:
            alle_list = list(map(float, alle_list))

        except:
            alle_list = list(map(str, alle_list))

        #  dictionary { allele : peak height }  ,   eg.{8.0: 499.0, 9.0: 105.0, 11.0: 246.0}
        dictionary = dict(zip(alle_list, height_list))
        alle_peakHeight_list.append(dictionary)

    #  dictionary { locus : { allele : peak height } }   eg.{D13S1358:{8.0: 499.0, 9.0: 105.0, 11.0: 246.0}}
    alle_peakHeight_dic = dict(zip(strloci_g, alle_peakHeight_list))

    return alle_peakHeight_dic


def adj_profile(tumor_profile_path):
    '''

    :param tumor_profile_path:  sample profile path
    :return:  the  adjusted  dictionary
    '''
    demo_profile = G_read.g_read(tumor_profile_path)
    logging.debug('profile dictionary: ')
    logging.debug(demo_profile)

    adj_profile = g_read_allePH(demo_profile)
    logging.debug('adjusted profile dictionary：')
    logging.debug(adj_profile)
    print('profile has been converted . ')
    return  adj_profile



