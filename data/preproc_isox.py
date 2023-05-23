# script for preprocessing data from orbitrap
# (c) by Anya for Ilya 22.May 2023

# run by typing in the command line inside the folder with files
#    # python3 preproc_isox.py 
# or in python notebook using 
#    # %run preproc_isox.py

# script reads all .isox files in the current folder as tables
# cuts tables into 4 separate ones
# resets time and scan number for each of the new tables
# writes down the files in the folder named by the original file

# if nesessary:
# adjust cutting times in the def main:
    # df1_split_times
    # df2_split_times
    # df3_split_times
    # df4_split_times

import glob
import os
import pandas as pd

def preproc_isox(path, df1_split_times,df2_split_times,df3_split_times,df4_split_times):

    #read in the file from the specified path
    data_isox = pd.read_table(path, sep='\s+')
    
    #print info
    print('The total length of the file in minutes is from \n' + str(data_isox['time.min'].min())+' min to '+
      str(data_isox['time.min'].max()) + ' min')
    print('____________________________________')
    print('The chosed split of the file is into 4 files \n 1st ' 
        + str(df1_split_times) +
        ' min \n 2nd ' + str(df2_split_times) +
        ' min \n 3rd ' + str(df3_split_times) +
        ' min \n 4th ' + str(df4_split_times) +' min')

    #split the dataframe into 4 based on provided data
    df1 = data_isox[(data_isox['time.min'] > df1_split_times[0]) & (data_isox['time.min'] < df1_split_times[1])] 
    df2 = data_isox[(data_isox['time.min'] > df2_split_times[0]) & (data_isox['time.min'] < df2_split_times[1])]
    df3 = data_isox[(data_isox['time.min'] > df3_split_times[0]) & (data_isox['time.min'] < df3_split_times[1])]
    df4 = data_isox[(data_isox['time.min'] > df4_split_times[0]) & (data_isox['time.min'] < df4_split_times[1])]
    
    # reset times and scan number in the new dataframes
    dfs_list=[df1,df2,df3,df4]
    for df in dfs_list:
        df.loc[:,'time.min']=df['time.min']-df['time.min'].min()+0.01 
        df.loc[:,'scan.no']=df['scan.no']-df['scan.no'].min()+1

    # #write down the dataframes
    print('____________________________________')
    print('Current file is ' + path)
    print('____________________________________')
    print('Currently saving:')
    
    for df in dfs_list:
        outpath=df['filename'].iloc[0]
        os.makedirs(outpath, exist_ok=True)  
        name=(outpath +'/'+ df['filename'].iloc[0]+'_'+str(df['compound'].min())+'_'+str(df['compound'].max())+'.isox')
        
        print(name)
        df.to_csv(name, sep='\t',float_format='%.3f')
    print('____________________________________')
    print('New files are saved in ./' + outpath +'/')

def main():    
    #path=glob.glob(os.path.join('../data/',"*.isox"))
    paths=glob.glob(os.path.join("*.isox"))

    # Set split times for the new files
    df1_split_times=[0,120]
    df2_split_times=[120.05,210]
    df3_split_times=[210.05,300]
    df4_split_times=[300.05,330]

    # these lines runs the function
    for path in paths:
        # do not process if folder exists 
        fold_name=path.split('.')[0]
        isExist=os.path.exists('./'+fold_name)
        if not isExist:
            print('does not exist -> start processing')
            preproc_isox(path, df1_split_times,df2_split_times,df3_split_times,df4_split_times)
        else:
            print(fold_name +' already exists')
            continue

if __name__ == "__main__":
    main() 