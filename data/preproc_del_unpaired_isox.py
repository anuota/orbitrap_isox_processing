# script for preprocessing data from orbitrap
#    # 1) splits the file based on times provided
#    # 2) deletes unpaired compounds for each scan
# (c) by Anya for Ilya 22.June 2023 with edits by Denis

# run by typing in the command line inside the folder with files
#    # python3 preproc_del_unpaired_isox.py 
# or in python notebook using 
#    # %run preproc_del_unpaired_isox.py 

# script reads all .isox files in the current folder as tables
# cuts tables into 4 separate ones
# resets time and scan number for each of the new tables
# writes down the files in the folder named by the original file

# if nesessary:
# adjust cutting times in the def main:
#   df1_split_times
#   df2_split_times
#   df3_split_times
#   df4_split_times

import glob
import os
import pandas as pd


def preproc_isox(path, df1_split_times, df2_split_times, df3_split_times, df4_split_times):
    # read in the file from the specified path
    data_isox = pd.read_table(path, sep='\s+')

    # print info
    print('The total length of the file is ' + str(data_isox['time.min'].min())+' min to '+
      str(data_isox['time.min'].max()) + ' min')
    print('____________________________________')
    print('The chosed split of the file is into 4 files \n 1st ' 
        + str(df1_split_times) +
        ' min \n 2nd ' + str(df2_split_times) +
        ' min \n 3rd ' + str(df3_split_times) +
        ' min \n 4th ' + str(df4_split_times) +' min')

    # split the dataframe into 4 based on provided data
    df1 = data_isox[(data_isox['time.min'] > df1_split_times[0]) & (data_isox['time.min'] < df1_split_times[1])] 
    df2 = data_isox[(data_isox['time.min'] > df2_split_times[0]) & (data_isox['time.min'] < df2_split_times[1])]
    df3 = data_isox[(data_isox['time.min'] > df3_split_times[0]) & (data_isox['time.min'] < df3_split_times[1])]
    df4 = data_isox[(data_isox['time.min'] > df4_split_times[0]) & (data_isox['time.min'] < df4_split_times[1])]

    # reset times and scan number in the new dataframes
    dfs_list = [df1, df2, df3, df4]
    for df in dfs_list:
        df.loc[:, 'time.min'] = df['time.min'] - df['time.min'].min() + 0.01
        df.loc[:, 'scan.no'] = df['scan.no'] - df['scan.no'].min() + 1

    return dfs_list


def percent(part, whole):
    try:
        return 100 * float(part) / float(whole)
    except ZeroDivisionError:
        return 0


def find_prob_comp(dfs_list,fold_name,proc):
    # a function to find a problematic compound
    # dfs_list=[df1,df2,df3,df4]
    i = 1
    for df in dfs_list:
        print('___Working on file #' + str(i))
        with open("./output.txt", "a") as f:
            print('___Working on file #' + str(i),file=f)
        # find problematic compunds
        compounds = df['compound'].unique()
        # print(compounds)
        # print compound if its quantity does not match to number of scans
        prob_comps = [comp for comp in compounds if df[(df['compound'] == comp)]['isotopolog'].count() / 2 != df['scan.no'].max()]
        # print(prob_comps)

        # for each problmatic compound print number of scans and number of compounds   
        for prob_co in prob_comps:
            numm_scans = (df['scan.no'].max())
            num_comp = df[df['compound'] == prob_co]['isotopolog'].count() / 2
            # calculate percent of missing data
            perc = 100 - percent(num_comp, numm_scans)
            if perc>=proc:
                df=df[~(df['compound'] == prob_co)]
            
            with open("./output.txt", "a") as f:
                if perc>=proc:
                    print('Missing more than ' +str(proc)+' data on compound '+str(prob_co)+'. Deleted ALL.',file=f)
                if perc>=3:
                    print('Problematic compound: ' + str(prob_co) + '; No of scans/compounds: ' + str(numm_scans) + ' / ' + str(num_comp), file=f)
                    print(f'   missing {perc:.2f} % of the data', file=f)
        i = i + 1

def delete_unpaired(dfs_list):
    # delete scans of compounds with unpaired isotopologs
    print('Starting to delete unpaired lines')
    dfs_list_new = []
    for df in dfs_list:
        # group by scan number
        pairs = df.groupby(['scan.no', 'compound'])['isotopolog'].size().reset_index(name='counts')

        # data to keep
        pairs_to_keep = pairs[pairs['counts'] == 2].drop(columns='counts')

        # keep only nesessary data in dataframe
        df = df.merge(pairs_to_keep[['scan.no', 'compound']], on=['scan.no', 'compound'], how='inner')

        dfs_list_new.append(df)
        print('Unpaired lines are deleted')
    return dfs_list_new


def save(dfs_list,fold_name):
    # # #write down the dataframes
    print('____________________________________')
    print('Currently saving:')
    i=1
    for df in dfs_list:
        #outpath=df['filename'].iloc[0]
        outpath=fold_name
        #print(outpath)
        os.makedirs(outpath, exist_ok=True)  
        #name=(outpath +'/'+ df['filename'].iloc[0]+'_'+str(df['compound'].min())+'_'+str(df['compound'].max())+'.isox')
        name=(outpath +'/'+ outpath+'_pt'+str(i)+'.isox')
        
        print(name)
        i=i+1
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

    proc=15 # min peccent of missing data to delete ALL info about compound

    # these lines runs the function
    for path in paths:
        # do not process if folder exists 
        fold_name=path.split('.')[0]
        isExist=os.path.exists('./'+fold_name)
        if not isExist:
            print(f'Folder \'{fold_name}\' does not exist -> start processing')
            #preprocessing
            dfs_list = preproc_isox(path, df1_split_times,df2_split_times,df3_split_times,df4_split_times)
            print('____________________________________')
            print('Current file is ' + path)
            with open("./output.txt", "a") as f:  
                print('____________________________________',file=f)
                print('Current file is ' + path,file=f)
            #if del #put a condition depending on the input
            find_prob_comp(dfs_list,fold_name,proc)
            #delete unpaired scans with unpaired isotopologs
            dfs_list = delete_unpaired(dfs_list)

            #save files
            save(dfs_list,fold_name)
        else:
            print(f'Folder \'{fold_name}\' already exists')
            continue

if __name__ == "__main__":
    main() 