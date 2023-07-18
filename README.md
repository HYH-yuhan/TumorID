# TumorID
TumorID, short for "tumor identification", is a supplementary for the manuscript entitled "Exploration of identifying individual tumor tissue using probabilistic genotyping algorithms". 

## TumorID_code
The pipeline has ability to calculate data containing four cell populations.  
Source code in this file is written by python 3.8.8.  **tumor_mian_iter.py** can be used as the mian program. Users can calculate the maximum likelihood under hp and hd, and store results in specified path by adjust the following paths : 
```
sample_name = '0166'
tumor_profile_path = './HTFD_%s.txt' %(sample_name)                   
sus_path = r'./HBD_%s.txt' %(sample_name)                             
fre_path = r'./southeast_fre.csv'                                         
sizepath = r'./gf_kit.txt' 
tmp_path = r'./example/'  
```
Each of them is the path to store tumor sample STR profile, reference profile, frequency table and size table of kits. Specially, tmp_path in this script is used to store each locus genotype combination file during the running process, which is created by users.  

**G_read.py** and **G_read_allePH.py** are used to convert the 'profile'.txt  into a dictionary; **tumor_geno_comb.py** is used to generate genotype combinations of each locus and put it in a temporary folder(tmp_path); **cal_P_g_H.py** is used to calculate  f(Hd) 、f(Hp) of Genotype combinations at specific locus; 
**P_g_H_dict.py** is used to put f(Hd) 、f(Hp) at different loci in a dictionary for ease of calculation; **cal_weighting.py** is used to calculate weighting , involving unknown parameters; **cal_h_product.py** is used to calculate likelihood values of all loci under hd / hp; **tumor_optimize.py** is used to search parameter and calculate likelihood; **tumor_main_iter.py** is for file input and results output. There are two result files for the likelihood values under the two hypotheses (Hp and Hd).  **info_read.py**,**initial_value.py**,**OSprocess.py** are  other function files that need to be used throughout the pipeline.  

For non-contributor testing operations, **Generate_p.py** and **Generate_pfile.py** are used. Firstly, designate allele list files generated using frequency files as script inputs(alle_path). Secondly, set the output path of the results ,which is a simulation-generated non-contributor genotyping file. Finally, use this file to replace the sus_path above and perform the calculations.
```
for i in range(1,1001):
    alle_path = r'./southeast_allele.txt'
    RM_path = r'/{}.txt'.format(i)
    Generate_p.RM_txt_simu(alle_path,RM_path)
```
  
For kinship consideration, list kls in the script **cal_P_g_H.py** should be adjusted.
```
 kls = [ 1 , 0 , 0 ]
```

Please note that the script is only available for globalfiler. Please contact the author of the manuscript if you encounter any bugs.

## Example
An example is provided here, including the **tumor sample profile HTFD0166.txt, reference profile HBD0166.txt, frequency file southeast_fre.csv, size file gf_kit.txt**. All of profile were exported by GeneMapper ID-X 1.5 software (Applied Biosystems, USA), frequency file were derived from a population survey of Southwest Han Chinese individuals, and globalfiler kit info was exported from [euroformix](http://euroformix.com/) **southeast_allele.txt** is used to generate random unrelated person.

## ValiOpt_dataset
This dataset is used to verify the accuracy of the optimizer in this pipeline compared with Euroformix. Eight two-person DNA mixture profiles, reference profile, and frequency file are provided. 

## Info about needed python-packages
numpy  : mathematical operation over arrays  
logging: generate files to record genotype combinations, genotype probabilities, and likelihood values for each locus  
scipy  : used for gamma distribution  
pandas : construct dataframe for imported evidence files, reference files, frequency tables, etc.  
os     : used for path Management
## Info about needed python-modules
itertools: functional tools for creating and using iterators.  
operator/functools : used for dimensionality reduction of lists 



* Erhii and Haoyu Wang also contributed to this script.
