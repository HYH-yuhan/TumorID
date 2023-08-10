# -*- coding: utf-8 -*-
"""
@author  : HYH
@E-mail  : huyuhan28@163.com
@className : P_g_H_dict
@describe  : TODO  put f(Hd) 、f(Hp) at different loci in a dictionary for ease of calculation
"""

import cal_P_g_H

def h_f_dict(out_path, fre_path, profile,sus_path,noc):
    '''

    :return: two dictionary，The former is the result of frequency calculation under hp,
    the latter is the result of frequency calculation under hd .
    eg {'D3S1358': [0.03580327524328875, 0.0036276920496056252, ...

    '''
    dict_fp = {}
    dict_fd = {}

    loci_ls = list(profile.keys())

    loci_ls.remove('AMEL')              ### Need to remove sexchromosome-related loci
    loci_ls.remove('DYS391')
    loci_ls.remove('Yindel')


    for locus in loci_ls:

        if locus == 'D3S1358':
            fai = 0.093
        if locus == 'vWA':
            fai = 0.1667
        if locus == 'D16S539':
            fai = 0.093
        if locus == 'CSF1PO':
            fai = 0.2326
        if locus == 'TPOX':
            fai = 0.0775
        if locus == 'D8S1179':
            fai = 0.1628
        if locus == 'D21S11':
            fai = 0.1473
        if locus == 'D18S51':
            fai = 0.4341
        if locus == 'D2S441':
            fai = 0.1685
        if locus == 'D19S433':
            fai = 0.1473
        if locus == 'TH01':
            fai = 0.1008
        if locus == 'FGA':
            fai = 0.2403
        if locus == 'D22S1045':
            fai = 0.1685
        if locus == 'D5S818':
            fai = 0.2481
        if locus == 'D13S317':
            fai = 0.1328
        if locus == 'D7S820':
            fai = 0.1008
        if locus == 'SE33':
            fai = 0.1685
        if locus == 'D10S1248':
            fai = 0.1685
        if locus == 'D1S1656':
            fai = 0.1685
        if locus == 'D12S391':
            fai = 0.1395
        if locus == 'D2S1338':
            fai = 0.1318

        g_fp_ls,g_fd_ls = cal_P_g_H.p_g_h(out_path,fre_path,profile,locus,fai,sus_path,noc)
        dict_fp[locus] = g_fp_ls
        dict_fd[locus] = g_fd_ls

    return dict_fp,dict_fd

