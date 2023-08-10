# -*- coding: utf-8 -*-
"""
@author  : Forenwhy / HYH
@E-mail  : wanghy0707@gmail.com / huyuhan28@163.com
@className : OSprocess
"""

import os

def mkResultFile (path):
    global result_file
    result_file = open(path, 'a')
    file_header = ['SampleName','Mx','Likelihood', 'Success', 'Reason', 'Mu', 'Omega','deg','Mx2','RunningTime']
    result_file.write("\t".join(str(b) for b in file_header))
    result_file.write('\n')

def listFilename(path):
    '''
    List all the file names under a folder
    '''
    filename_list = []
    for info in os.listdir(path):
        filename_list.append(info)
    return filename_list

def listFilename_nosuffix(path):
    '''
    List all the file names under a folder (excluding .txt)
    '''
    listFilename_nosuffix_ls =[]
    for info in os.listdir(path):
        info = info[0:-4]
        listFilename_nosuffix_ls.append(info)
    return listFilename_nosuffix_ls


def openFile(path, separate_sym):
    listforfile = []
    if separate_sym != None:
        with open(path, 'r') as file:
            for line in file:
                listforfile.append(list(line.strip('\n').split(separate_sym)))
    elif separate_sym == None:
        with open(path, 'r') as file:
            for line in file:
                listforfile.append(line.strip('\n'))
    return listforfile

