# TumorID
TumorID, short for "tumor identification", is a supplementary for the manuscript entitled "Exploration of individual identification of tumor tissues using probabilistic genotyping algorithms". 

## TumorID_code
Source code in this file is written by python 3.8.8. tumor_mian_iter.py can be used as the mian program. Users can calculate the maximum likelihood under hp and hd, and store results in specified path by adjust the following paths : 

```
sample_name = '0166'
tumor_profile_path = './HTFD_%s.txt' %(sample_name)                   
sus_path = r'./HBD_%s.txt' %(sample_name)                             
fre_path = r'./southeast_fre.csv'                                         
sizepath = r'./gf_kit.txt' 
tmp_path = r'./example/'  
```
Each of them are the path to store tumor sample STR profile, reference profile, frequency table and size table of kits. Specially, tmp_path in this script is used to store each locus genotype combination file during the running process, which is created by users.  

Please note that the script is only available for globalfiler.Please contact the author of the manuscript if you encounter any bugs.
## Example
An example is provided here, including the **tumor sample profile HTFD0166.txt, reference profile HBD0166.txt, frequency file southeast_fre.txt, size file gf_kit.txt**. All of profile were exported by GeneMapper ID-X 1.5 software (Applied Biosystems, USA), frequency file were derived from a population survey of Southwest Han Chinese individuals,and globalfiler kit info was exported from [euroformix](http://euroformix.com/)

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
