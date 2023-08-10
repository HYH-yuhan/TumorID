# -*- coding: utf-8 -*-
"""
@author  : HYH
@E-mail  : huyuhan28@163.com
@className : cal_weighting
@describe  : TODO calculate weighting ,involving unknown parameters
"""


import numpy as np
from scipy.stats import gamma
import cal_P_g_H
from info_read import read_size

def genocom_numerator_ls(mx,x,out_path,profile,locus,sizefile,noc):
    '''

    :param mx:  mixture proportion of the normal cell.
    :param x:  unknown parameters, noc = 2  x = [mu , omega , deg] ；
                                    noc = 3  x = [mu , omega , deg , proportion of the second tumor cell population]；
                                    noc = 4  x = [mu , omega , deg , proportion of the second / third tumor cell population ]
    :param sizefile: size dataframe
    :return: Numerator of the weights under each genotype combination,considering mx and degradation parameter,Two-dimensional list
            [[1.8, 0.19999999999999996, 0]
    '''

    ls = cal_P_g_H.get_genocom_ls(out_path, locus)

    locus_allele_q  = cal_P_g_H.get_allele_ls_Q(profile,locus)

    genocom_numerator_ls = []
    if noc == 2:
        for info in ls:  # for each genotype combination

            geno_numeratordg_ls = []
            geno_numerator_ls = []  # [1.8, 0.19999999999999996, 0]

            i = 0
            while i < len(locus_allele_q):
                numerator = 0

                if info[0] == locus_allele_q[i]:
                    numerator += mx

                if info[1] == locus_allele_q[i]:
                    numerator += mx

                if info[2] == locus_allele_q[i]:
                    numerator = numerator + 1 - mx

                if info[3] == locus_allele_q[i]:
                    numerator = numerator + 1 - mx

                geno_numerator_ls.append(numerator)
                i += 1

            i = 0
            for numera in geno_numerator_ls:
                Size = read_size(sizefile, locus, locus_allele_q[i])
                adjSize = (Size - 90) / 100
                numera = numera * x[2] ** adjSize
                geno_numeratordg_ls.append(numera)
                i += 1
            genocom_numerator_ls.append(geno_numeratordg_ls)

    # Numerator
    if noc == 3:
        for info in ls:

            geno_numeratordg_ls = []

            geno_numerator_ls = []

            i = 0

            while i < len(locus_allele_q):
                numerator = 0
                if info[0] == locus_allele_q[i]:
                    numerator += mx

                if info[1] == locus_allele_q[i]:
                    numerator += mx

                if info[2] == locus_allele_q[i]:
                    numerator = numerator + x[3]

                if info[3] == locus_allele_q[i]:
                    numerator = numerator + x[3]

                if info[4] == locus_allele_q[i]:
                    numerator = numerator + (1 - mx - x[3])

                if info[5] == locus_allele_q[i]:
                    numerator = numerator + (1 - mx - x[3])

                geno_numerator_ls.append(numerator)
                i += 1

            i = 0
            for numera in geno_numerator_ls:
                Size = read_size(sizefile, locus, locus_allele_q[i])
                adjSize = (Size - 90) / 100
                numera = numera * x[2] ** adjSize

                geno_numeratordg_ls.append(numera)
                i += 1
            genocom_numerator_ls.append(geno_numeratordg_ls)

    if noc == 4:
        for info in ls:

            geno_numeratordg_ls = []
            geno_numerator_ls = []  # [1.8, 0.19999999999999996, 0]

            i = 0

            while i < len(locus_allele_q):
                numerator = 0
                if info[0] == locus_allele_q[i]:
                    numerator += mx

                if info[1] == locus_allele_q[i]:
                    numerator += mx

                if info[2] == locus_allele_q[i]:
                    numerator = numerator + x[3]

                if info[3] == locus_allele_q[i]:
                    numerator = numerator + x[3]

                if info[4] == locus_allele_q[i]:
                    numerator = numerator + x[4]

                if info[5] == locus_allele_q[i]:
                    numerator = numerator + x[4]

                if info[6] == locus_allele_q[i]:
                    numerator = numerator + (1 - mx - x[3] - x[4])

                if info[7] == locus_allele_q[i]:
                    numerator = numerator + (1 - mx - x[3] - x[4])

                geno_numerator_ls.append(numerator)
                i += 1

            i = 0
            for numera in geno_numerator_ls:
                Size = read_size(sizefile, locus, locus_allele_q[i])
                adjSize = (Size - 90) / 100
                numera = numera * x[2] ** adjSize
                geno_numeratordg_ls.append(numera)
                i += 1

            genocom_numerator_ls.append(geno_numeratordg_ls)
    return  genocom_numerator_ls

def genocom_weighting_ls(mx,x,AT,out_path,profile,locus,sizefile,noc):
    '''

        :return: Weighting of each allele under each genotype combination, two-dimensional list
                [[2.57490496e-04 1.39455040e-04 1.00000000e+00]。。。

    '''
    locus_allele_q = cal_P_g_H.get_allele_ls_Q(profile,locus)

    genocom_weighting_ls = []

    for info in genocom_numerator_ls(mx,x,out_path,profile,locus,sizefile,noc):
        alpha_ls = [i/(x[1]*x[1]) for i in info ]

        weighting_ls = []
        for alle, alpha in zip(locus_allele_q, alpha_ls):
            if alle != 'Q':
                weighting = gamma.pdf(profile[locus][float(alle)], alpha, scale=(x[0] * x[1] * x[1]))

            else:
                weighting = gamma.cdf(AT, alpha, scale=x[0] * x[1] * x[1])
            weighting_ls.append(weighting)
        genocom_weighting_ls.append(weighting_ls)


    genocom_weighting_ls = np.nan_to_num(genocom_weighting_ls, nan=1)
    #logging.debug('weighting list of each allele :  ')
    #logging.debug(genocom_weighting_ls)
    return genocom_weighting_ls

def weighting(mx,x,AT,out_path,profile,locus,sizefile,noc):
    '''

    :return: Weighting product of all genotype combinations, list

    '''
    weighting_ls = []

    for geno in genocom_weighting_ls(mx,x,AT,out_path,profile,locus,sizefile,noc):
        result = np.prod(geno)
        weighting_ls.append(result)
    #logging.debug('weighting list of genotype combination:  ')
    #logging.debug(weighting_ls)
    return weighting_ls