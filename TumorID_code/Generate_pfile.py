# -*- coding: utf-8 -*-
"""
@author  : HYH
@E-mail  : huyuhan28@163.com
@date    : 2023 / 04 / 05  19:37:26
@className : Generate_pfile
@describe  : TODO unrelated individual simulation
"""
import  Generate_p

for i in range(1,1001):
    alle_path = r'.//southeast_allele.txt'
    RM_path = './/{}.txt'.format(i)
    Generate_p.RM_txt_simu(alle_path,RM_path)