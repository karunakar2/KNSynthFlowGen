import numpy as np
import pandas as pd

import sys
import os

import glob

#from base import *
import base

#for parallel execution, this file has to be named main.py
#else it will break the parallel part of the code
#https://docs.python.org/3/library/multiprocessing.html#multiprocessing-programming
name = __name__ 

def main():
    DaysPerMonth = [31,28,31,30,31,30,31,31,30,31,30,31]    
    
    qDaily = pd.read_csv('../data/Qdaily_.txt', parse_dates=["Date"]) 
    #this file has dateinfo
    qDaily.set_index('Date',inplace=True)
    for thisStn in qDaily.columns:
        if 'evap' in thisStn:
            qDaily[thisStn] = np.exp(qDaily[thisStn])

    nYears = int(qDaily.index.year.max() - qDaily.index.year.min())+1
    nSites = len(qDaily.columns)
    #print(qDaily.info())

    #my specs
    ##nR = min(100, pow(nYears,2))
    ##num_realizations = [nYears, nR] #lall's criterion
    ##num_years = [nYears, 1]
    
    #matlab specs
    num_realizations = [100, 1000]
    num_years = [100, 1]
    
    #test set
    #num_realizations = [3]
    #num_years = [2]
    
    #set folders
    try:
        os.makedirs('../validation', exist_ok=True)
    except:
        pass
        
    try:
        os.makedirs('../validation/synthetic/', exist_ok=True)
        #os.chdir('../validation/synthetic/')
    except Exception as err:
        raise Exception(err)
    
    for fName in glob.glob('../validation/synthetic/*.csv'):
        try:
            os.remove(fName)
        except Exception as er:
            print(er)
            pass
    
    #Kirsch + Nowak generation
    for r,y in zip(num_realizations,num_years): 
        print(r,'realisations -',y,'years')
        Qd_cg, _ = base.combined_generator(qDaily, r, y, name=name)

        for i in Qd_cg.keys(): #realisations
            #daily
            for thisStn in Qd_cg[i].columns: #stations
                if 'evap' in thisStn:
                    # back-transform evaporation
                    Qd_cg[i][thisStn] = np.log(Qd_cg[i][thisStn].values)
                fName = '../validation/synthetic/' + str(thisStn)+'-'+str(r)+'x'+str(y) + '-daily.csv'
                with open(fName,'a') as f:
                    f.write(','.join(Qd_cg[i][thisStn].astype('str'))+ '\n')
            
            #compounded file
            Qd_cg[i].to_csv('./../validation/synthetic/Qdaily-'+str(r)+'x'+str(y)+'.csv', mode='a', header=False, index=False)
                        
            #monthly
            monDf = base.convert_data_to_monthly(Qd_cg[i]) #, debug=True)
            for thisStn in monDf.keys():
                if 'evap' in thisStn:
                    # divide evaporation by 86400 (s/day) to get total monthly evap in mm/month
                    #monDf[thisStn] = monDf[thisStn].divide(pd.Series(DaysPerMonth),index=monDf.columns) #daily conversion
                    monDf[thisStn] /= 86400 #this is an artifact in the monthly fn
                    #print(thisStn)
                temp = np.ravel(monDf[thisStn].to_numpy())
                fName = '../validation/synthetic/' + str(thisStn)+'-'+str(r)+'x'+str(y) + '-monthly.csv'
                with open(fName,'a') as f:
                    temp.tofile(f, sep=',')
                    f.write('\n')
            #monDf[thisStn].to_csv('../validation/synthetic/r'+str(i)+'X-y'+str(y)+thisStn+'-monthly.csv')

if __name__ == "__main__":
    main()






