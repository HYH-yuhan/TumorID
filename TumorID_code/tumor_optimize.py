# -*- coding: utf-8 -*-
"""
@author  : HYH
@E-mail  : huyuhan28@163.com
@className : tumor_optimize
@describe  : TODO  The main program of the parameter searching
"""

import cal_h_product
from scipy.optimize import minimize   # Mention: Scipy >1.7.1
import logging
import time
import initial_value
import numpy as np
import warnings


#  the likelihood formula under hp ,x is variable
def fun_hp(mx, AT,tmp_path, sus_path, profile,dict_fp,sizefile,noc):
    v = lambda x:  cal_h_product.products_hp(mx, x, AT,tmp_path, sus_path, profile,dict_fp,sizefile,noc)
    logging.debug('-----------------------------during the likelihood calculation----------------------------------------')
    return v

# the likelihood formula under hd ,x is variable
def fun_hd(mx, AT,tmp_path, profile,dict_fd,sizefile,noc):
    v = lambda x:  cal_h_product. products_hd(mx, x, AT,tmp_path, profile,dict_fd,sizefile,noc)
    logging.debug('-----------------------------during the likelihood calculation----------------------------------------')
    return v

# limitation range of unknown parameter
def get_cons(noc,mx):
    if noc == 2:
        args2 =  (0.1,20000.0 , 0.00001,1.0 , 0.00001,1.0)
        x0min, x0max, x1min, x1max, x2min, x2max = args2
        cons = ({'type': 'ineq', 'fun': lambda x: x[0] - x0min},
                {'type': 'ineq', 'fun': lambda x: -x[0] + x0max},
                {'type': 'ineq', 'fun': lambda x: x[1] - x1min},
                {'type': 'ineq', 'fun': lambda x: -x[1] + x1max} ,
                {'type': 'ineq', 'fun': lambda x: x[2] - x2min},
                {'type': 'ineq', 'fun': lambda x: -x[2] + x2max})
    elif noc == 3 :
        args3 = (0.1, 20000.0, 0.00001, 1.0, 0.00001, 1.0, 0.00001, 1 - mx)
        x0min, x0max, x1min, x1max, x2min, x2max, x3min, x3max = args3
        cons = ({'type': 'ineq', 'fun': lambda x: x[0] - x0min},
                {'type': 'ineq', 'fun': lambda x: -x[0] + x0max},
                {'type': 'ineq', 'fun': lambda x: x[1] - x1min},
                {'type': 'ineq', 'fun': lambda x: -x[1] + x1max},
                {'type': 'ineq', 'fun': lambda x: x[2] - x2min},
                {'type': 'ineq', 'fun': lambda x: -x[2] + x2max},
                {'type': 'ineq', 'fun': lambda x: x[3] - x3min},
                {'type': 'ineq', 'fun': lambda x: -x[3] + x3max})

    elif noc == 4 :
        args4 = (0.1, 20000.0, 0.00001, 1.0, 0.00001, 1.0, 0.0000001, 1 - mx)
        x0min, x0max, x1min, x1max, x2min, x2max, xmin, xmax = args4
        cons = ({'type': 'ineq', 'fun': lambda x: x[0] - x0min},
                {'type': 'ineq', 'fun': lambda x: -x[0] + x0max},
                {'type': 'ineq', 'fun': lambda x: x[1] - x1min},
                {'type': 'ineq', 'fun': lambda x: -x[1] + x1max},
                {'type': 'ineq', 'fun': lambda x: x[2] - x2min},
                {'type': 'ineq', 'fun': lambda x: -x[2] + x2max},
                {'type': 'ineq', 'fun': lambda x: x[3] - xmin},
                {'type': 'ineq', 'fun': lambda x: x[4] - xmin},
                {'type': 'ineq', 'fun': lambda x: -x[3] + xmax},
                {'type': 'ineq', 'fun': lambda x: -x[4] + xmax},
                {'type': 'ineq', 'fun': lambda x: -x[3]-x[4] + xmax})
    return cons

# obtain initial Value
def get_start_val (noc,tumor_profile_path):
    start_val = initial_value.ini_val(tumor_profile_path)
    if noc == 2:
        x0 = np.asarray((start_val[0], start_val[1], 0.5))
    elif noc == 3:
        x0 = np.asarray((start_val[0], start_val[1], 0.5, 0.1))
    elif noc == 4:
        x0 = np.asarray((start_val[0], start_val[1], 0.5, 0.1, 0.1))
    return x0


# likelihood calculation and result output
def optFuction(fun,tumor_profile_path,sample_name,mx,outt_path,noc):
    start_time = time.time()
    x0 = get_start_val (noc,tumor_profile_path)
    cons = get_cons(noc,mx)

    logging.debug('--------------------------------Start likelihood calculation-----------------------------------------')
    warnings.filterwarnings("ignore", category=RuntimeWarning)

    try:
        res = minimize(fun, x0, method='Nelder-Mead', constraints = cons)
    except KeyError as e:
        print(e)
        return

    logging.debug('+1')
    end_time = time.time()
    cost_time = end_time - start_time
    if len(x0) == 3:
        re_ls = [sample_name, mx, res.fun, res.success, res.message, res.x[0], res.x[1], res.x[2], cost_time]
    elif len(x0) == 4 :
        re_ls = [sample_name, mx, res.fun, res.success, res.message, res.x[0], res.x[1], res.x[2],res.x[3] ,cost_time]
    elif len(x0) == 5 :
        re_ls = [sample_name, mx, res.fun, res.success, res.message, res.x[0], res.x[1], res.x[2], res.x[3],  res.x[4],cost_time]

    result_file = open(outt_path, 'a')
    result_file.write("\t".join(str(b) for b in re_ls))
    result_file.write('\n')






