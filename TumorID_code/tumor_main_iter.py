# -*- coding: utf-8 -*-
"""
@author  : HYH
@E-mail  : huyuhan28@163.com
@className : tumor_main_iter
@describe  : TODO file input and result output
"""
import pandas as pd
import P_g_H_dict
import G_read_allePH
import tumor_geno_comb
import logging
import numpy as np
import initial_value
from OSprocess import mkResultFile
from tumor_optimize import fun_hp,fun_hd,optFuction

# sample
sample_name = '0166'
logging.basicConfig(filename='myprogramLog%s.txt'%(sample_name), level=logging.DEBUG,
                    format=' %(asctime)s - %(levelname)s- %(message)s')

# path
tumor_profile_path = './/HTFD_%s.txt' %(sample_name)                   #   path , sample profile
sus_path = r'.//HBD_%s.txt' %(sample_name)                             #   path , reference profile
fre_path = r'.//southeast_fre.csv'                                          #   path , frequency table
sizepath = r'.//gf_kit.txt'                                                 #   path , size table
tmp_path = r'.//example/'                                                     #   path to store each locus genotype combination file
outt_path_Hp = './/HTFD%s_result_Hp.txt' % (sample_name)                    #   path , likelihood result file under hp
outt_path_Hd = './/HTFD%s_result_Hd.txt' % (sample_name)                    #   path , likelihood result file under hd

sizefile = pd.read_table(sizepath)
#  Creating the result file
mkResultFile(outt_path_Hp)
mkResultFile(outt_path_Hd)

#   priori parameters
mx_ls = np.linspace(0.3,0.5,21)  # or mx_ls = [0.4]
AT = 175.0
noc = initial_value.noc_val(tumor_profile_path)


for mx in mx_ls:

    #  【first】 Generate genotype combination files for different loci based on the sample profile
    profile = G_read_allePH.adj_profile(tumor_profile_path)
    logging.debug(profile)
    tumor_geno_comb.profile_geno_comb(noc, tmp_path, profile)

    #  【second】 Generate frequency dictionary
    dict_fp, dict_fd = P_g_H_dict.h_f_dict(tmp_path, fre_path, profile,sus_path,noc)
    logging.debug(dict_fp)
    logging.debug(dict_fd)

    #   x: unknown parameters, noc = 2  x = [mu , omega , deg] ；
    #                          noc = 3  x = [mu , omega , deg , proportion of the second tumor cell population]；
    #                          noc = 4  x = [mu , omega , deg , proportion of the second / third tumor cell population ]
    #  【third】  likelihood calculation
    optFuction(fun_hp(mx, AT,tmp_path, sus_path, profile,dict_fp,sizefile,noc), tumor_profile_path,  sample_name, mx,
               outt_path_Hp,noc)
    optFuction(fun_hd(mx, AT,tmp_path, profile,dict_fd,sizefile,noc), tumor_profile_path, sample_name, mx, outt_path_Hd,noc)

    ##  when conduct non-contributor test, if ref cannot be interpreted without 'QQ', the fun_hp procedure is launched
    ##  there is no result in hp result profile , and the ""logging.info('Quit the program because the loci do not match： {}'.format(locus))"" can be founded in logging
