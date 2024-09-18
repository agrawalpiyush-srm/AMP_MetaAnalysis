######################################################################################
# AMPCell is developed for predicting, desigining and scanning anti-bacterial,       #
# anti-viral, and anti-fungal peptides.                                              #
######################################################################################
import argparse
import warnings
import subprocess
import textwrap
import os
import sys
import numpy as np
import pandas as pd
import math
import itertools
from collections import Counter
import pickle
import re
import glob
import time
import uuid
from time import sleep
from tqdm import tqdm
from sklearn.ensemble import ExtraTreesClassifier
import zipfile

warnings.filterwarnings('ignore', category=DeprecationWarning)

parser = argparse.ArgumentParser(description='Please provide following arguments', formatter_class=argparse.RawTextHelpFormatter) 

## Read Arguments from command
parser.add_argument("-i", "--input", type=str, required=True, help="Input: protein or peptide sequence(s) in FASTA format or single sequence per line in single letter code")
parser.add_argument("-o", "--output",type=str, help="Output: File for saving results by default outfile.csv")
parser.add_argument("-j", "--job",type=int, choices = [1,2,3], help="Job Type: 1:Predict, 2: Design, 3:Scan, by default 1")
parser.add_argument("-m", "--method",type=int, choices = [1,2,3], help="Model: 1:Anti-bacterial, 2: Anti-viral, 3:Anti-fungal, by default 1")
parser.add_argument("-t","--threshold", type=float, help=textwrap.dedent("""\Threshold: Value between 0 to 1 by default 0.32\nPlease use the following threshold values for best performance:\n1. For Anti-bacterial model: 0.32\n2. For Anti-viral model: 0.71\n3. For Anti-fungal model: 0.19"""))
parser.add_argument("-w","--winleng", type=int, choices =range(4, 101), help="Window Length: 4 to 100 (scan mode only), by default 4")
parser.add_argument("-d","--display", type=int, choices = [1,2], help="Display: 1:Positive Results Only, 2: All peptides, by default 1")
args = parser.parse_args()

# Function for generating all possible mutants
def mutants(file1,file2):
    std = list("ACDEFGHIKLMNPQRSTVWY")
    cc = []
    dd = []
    ee = []
    df2 = file2
    df2.columns = ['Name']
    df1 = file1
    df1.columns = ['Seq']
    for k in range(len(df1)):
        cc.append(df1['Seq'][k])
        dd.append('Original_'+'Seq'+str(k+1))
        ee.append(df2['Name'][k])
        for i in range(0,len(df1['Seq'][k])):
            for j in std:
                if df1['Seq'][k][i]!=j:
                    dd.append('Mutant_'+df1['Seq'][k][i]+str(i+1)+j+'_Seq'+str(k+1))
                    cc.append(df1['Seq'][k][:i] + j + df1['Seq'][k][i + 1:])
                    ee.append(df2['Name'][k])
    xx = pd.concat([pd.DataFrame(ee),pd.DataFrame(dd),pd.DataFrame(cc)],axis=1)
    xx.columns = ['Seq_ID','Mutant_ID','Seq']
    return xx
# Function for generating pattern of a given length
def seq_pattern(file1,file2,num):
    df1 = file1
    df1.columns = ['Seq']
    df2 = file2
    df2.columns = ['Name']
    cc = []
    dd = []
    ee = []
    for i in range(len(df1)):
        for j in range(len(df1['Seq'][i])):
            xx = df1['Seq'][i][j:j+num]
            if len(xx) == num:
                cc.append(df2['Name'][i])
                dd.append('Pattern_'+str(j+1)+'_Seq'+str(i+1))
                ee.append(xx)
    df3 = pd.concat([pd.DataFrame(cc),pd.DataFrame(dd),pd.DataFrame(ee)],axis=1)
    df3.columns= ['Seq_ID','Pattern_ID','Seq']
    return df3
# Function to check the seqeunce
def readseq(file):
    with open(file) as f:
        records = f.read()
    records = records.split('>')[1:]
    seqid = []
    seq = []
    for fasta in records:
        array = fasta.split('\n')
        name, sequence = array[0].split()[0], re.sub('[^ACDEFGHIKLMNPQRSTVWY-]', '', ''.join(array[1:]).upper())
        seqid.append('>'+name)
        seq.append(sequence)
    if len(seqid) == 0:
        f=open(file,"r")
        data1 = f.readlines()
        for each in data1:
            seq.append(each.replace('\n',''))
        for i in range (1,len(seq)+1):
            seqid.append(">Seq_"+str(i))
    df1 = pd.DataFrame(seqid)
    df2 = pd.DataFrame(seq)
    return df1,df2
# Function to check the length of seqeunces
def lenchk(file1):
    cc = []
    df1 = file1
    df1.columns = ['seq']
    for i in range(len(df1)):
        if len(df1['seq'][i])>100:
            cc.append(df1['seq'][i][0:100])
        else:
            cc.append(df1['seq'][i])
    df2 = pd.DataFrame(cc)
    df2.columns = ['Seq']
    return df2
# Function to generate the features out of seqeunces
def feature_gen(file,q=1):
    std = list('ACDEFGHIKLMNPQRSTVWY')
    df1 = file
    df1.columns = ['Seq']
    dd = []
    for j in df1['Seq']:
        cc = []
        for i in std:
            count = 0
            for k in j:
                temp1 = k
                if temp1 == i:
                    count += 1
                composition = (count/len(j))*100
            cc.append(composition)
        dd.append(cc)
    df2 = pd.DataFrame(dd)
    head = []
    for mm in std:
        head.append('AAC_'+mm)
    df2.columns = head
    return df2
def model_run(file1,file2):
    a = []
    data_test = file1
    clf = pickle.load(open(file2,'rb'))
    y_p_score1=clf.predict_proba(data_test)
    y_p_s1=y_p_score1.tolist()
    a.extend(y_p_s1)
    df = pd.DataFrame(a)
    df1 = df.iloc[:,-1].round(2)
    df2 = pd.DataFrame(df1)
    df2.columns = ['ML_score']
    return df2
def result_processor(name1,meth,ml_results,thresh):
        df2 = name1
        df2.columns = [0]
        df3 = ml_results
        ss = []
        for j in df2[0]:
            ss.append(j.replace('>',''))
        df6 = pd.DataFrame()
        df6['Seq_ID'] = ss
        df6['ML_Score'] = df3['ML_score']
        if meth ==1:
            df6['Prediction'] = ['Anti-bacterial' if df6['ML_Score'][i]>thresh else 'Non-anti-bacterial' for i in range(0,len(df6))]
        if meth ==2:
            df6['Prediction'] = ['Anti-viral' if df6['ML_Score'][i]>thresh else 'Non-anti-viral' for i in range(0,len(df6))]
        if meth ==3:
            df6['Prediction'] = ['Anti-fungal' if df6['ML_Score'][i]>thresh else 'Non-anti-fungal' for i in range(0,len(df6))]
        return df6    
print('############################################################################################')
######################################################################################
# AMPCell is developed for predicting, desigining and scanning anti-bacterial,       #
# anti-viral, and anti-fungal peptides.                                              #
######################################################################################
print('############################################################################################')

# Parameter initialization or assigning variable for command level arguments

Sequence= args.input        # Input variable 
 
# Output file 
 
if args.output == None:
    result_filename= "outfile.csv" 
else:
    result_filename = args.output
         
# Threshold 
if args.threshold == None:
        Threshold = 0.32
else:
        Threshold= float(args.threshold)
# Job Type 
if args.job == None:
        Job = int(1)
else:
        Job = int(args.job)
# Method Type
if args.method == None:
        met = int(1)
else:
        met = int(args.method)
# Window Length 
if args.winleng == None:
        Win_len = int(4)
else:
        Win_len = int(args.winleng)

# Display
if args.display == None:
        dplay = int(1)
else:
        dplay = int(args.display)


###########################################################################################

if Job==2:
    print("\n");
    print('#####################################################################################')
    print('Summary of Parameters:')
    print('Input File: ',Sequence,'; Threshold: ', Threshold,'; Job Type: ',Job,'; Method: ',met)
    print('Output File: ',result_filename,'; Window Length: ',Win_len,'; Display: ',dplay)
    print('#####################################################################################')
else:
    print("\n");
    print('######################################################################################')
    print('Summary of Parameters:')
    print('Input File: ',Sequence,'; Threshold: ', Threshold,'; Job Type: ',Job,'; Method: ',met)
    print('Output File: ',result_filename,'; Display: ',dplay)
    print('# ####################################################################################')
#======================= Prediction Module start from here =====================
if Job == 1:
    print('\n======= Thanks for using Predict module of AMP Cell. Your results will be stored in file :',result_filename,' =====\n')
    df_2,dfseq = readseq(Sequence)
    df1 = lenchk(dfseq)
    X = feature_gen(df1)
    if met == 1:
        mlres = model_run(X,'model/bac_model.sav')
    elif met == 2:
        mlres = model_run(X,'model/vir_model.sav')
    elif met == 3:
        mlres = model_run(X,'model/fun_model.sav')
    filename = str(uuid.uuid4())
    df11 = pd.concat([df_2,df1],axis=1)
    df11.to_csv(filename,index=None,header=False,sep="\n")
    mlres = mlres.round(3)
    df44 = result_processor(df_2,met,mlres,Threshold)
    df44['Sequence'] = df1.Seq
    df44 = df44[['Seq_ID','Sequence','ML_Score','Prediction']]
    if dplay == 1 and met==1:
        df44 = df44.loc[df44.Prediction=="Anti-bacterial"]
    elif dplay == 1 and met==2:
        df44 = df44.loc[df44.Prediction=="Anti-viral"]
    elif dplay == 1 and met==3:
        df44 = df44.loc[df44.Prediction=="Anti-fungal"]
    else:
        df44 = df44
    df44 = round(df44,3)
    df44.to_csv(result_filename, index=None)
    os.remove(filename)
    print("\n=========Process Completed. Have an awesome day ahead.=============\n")    
#===================== Design Model Start from Here ======================
elif Job == 2:
    print('\n======= Thanks for using Design module of AMP Cell. Your results will be stored in file :',result_filename,' =====\n')
    print('==== Designing Peptides: Processing sequences please wait ...')
    df_2,dfseq = readseq(Sequence)
    df1 = lenchk(dfseq)
    df_1 = mutants(df1,df_2)
    dfseq = df_1[['Seq']]
    X = feature_gen(dfseq)
    if met == 1:
        mlres = model_run(X,'model/bac_model.sav')
    elif met == 2:
        mlres = model_run(X,'model/vir_model.sav')
    elif met == 3:
        mlres = model_run(X,'model/fun_model.sav')
    filename = str(uuid.uuid4())
    df_1['Mutant'] = ['>'+df_1['Mutant_ID'][i] for i in range(len(df_1))]
    df11 = df_1[['Mutant','Seq']] 
    df11.to_csv(filename,index=None,header=False,sep="\n")
    mlres = mlres.round(3)
    df44 = result_processor(df11[['Mutant']],met,mlres,Threshold)
    df44['Mutant_ID'] = ['_'.join(df44['Seq_ID'][i].split('_')[:-1]) for i in range(len(df44))]
    df44.drop(columns=['Seq_ID'],inplace=True)
    df44['Seq_ID'] = [i.replace('>','') for i in df_1['Seq_ID']]
    df44['Sequence'] = df_1.Seq
    df44 = df44[['Seq_ID','Mutant_ID','Sequence','ML_Score','Prediction']]
    if dplay == 1 and met==1:
        df44 = df44.loc[df44.Prediction=="Anti-bacterial"]
    elif dplay == 1 and met==2:
        df44 = df44.loc[df44.Prediction=="Anti-viral"]
    elif dplay == 1 and met==3:
        df44 = df44.loc[df44.Prediction=="Anti-fungal"]
    else:
        df44 = df44
    df44 = round(df44,3)
    df44.to_csv(result_filename, index=None)
    os.remove(filename)
    print("\n=========Process Completed. Have an awesome day ahead.=============\n")
#=============== Scan Model start from here ==================
elif Job==3:
    print('\n======= Thanks for using Scan module of AMP Cell. Your results will be stored in file :',result_filename,' =====\n')
    print('==== Scanning Peptides: Processing sequences please wait ...')
    df_2,dfseq = readseq(Sequence)
    df_1 = seq_pattern(dfseq,df_2,Win_len)
    dfseq = df_1[['Seq']]
    X = feature_gen(dfseq)
    if met == 1:
        mlres = model_run(X,'model/bac_model.sav')
    elif met == 2:
        mlres = model_run(X,'model/vir_model.sav')
    elif met == 3:
        mlres = model_run(X,'model/fun_model.sav')
    filename = str(uuid.uuid4())
    df_1['Pattern'] = ['>'+df_1['Pattern_ID'][i] for i in range(len(df_1))]
    df11 = df_1[['Pattern','Seq']]
    df11.to_csv(filename,index=None,header=False,sep="\n")
    mlres = mlres.round(3)
    df44 = result_processor(df11[['Pattern']],met,mlres,Threshold)
    df44['Pattern_ID'] = ['_'.join(df44['Seq_ID'][i].split('_')[:-1]) for i in range(len(df44))]
    df44.drop(columns=['Seq_ID'],inplace=True)
    df44['Seq_ID'] = [i.replace('>','') for i in df_1['Seq_ID']]
    df44['Sequence'] = df_1.Seq
    df44 = df44[['Seq_ID','Pattern_ID','Sequence','ML_Score','Prediction']]
    if dplay == 1 and met==1:
        df44 = df44.loc[df44.Prediction=="Anti-bacterial"]
    elif dplay == 1 and met==2:
        df44 = df44.loc[df44.Prediction=="Anti-viral"]
    elif dplay == 1 and met==3:
        df44 = df44.loc[df44.Prediction=="Anti-fungal"]
    else:
        df44 = df44
    df44 = round(df44,3)
    df44.to_csv(result_filename, index=None)
    os.remove(filename)
    print("\n=========Process Completed. Have an awesome day ahead.=============\n")
print('\n======= Thanks for using AMP Cell. Your results are stored in file :',result_filename,' =====\n\n')
