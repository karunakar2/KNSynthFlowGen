import numpy as np
import pandas as pd

#from base import *
import base

#for parallel execution, this file has to be named main.py
#else it will break the parallel part of the code
#https://docs.python.org/3/library/multiprocessing.html#multiprocessing-programming
name = __name__ 

def main():
    qDaily = pd.read_csv('../data/Qdaily.txt', parse_dates=["Date"])
    qDaily.set_index('Date',inplace=True)
    for thisStn in qDaily.columns:
        if 'evap' in thisStn:
            qDaily[thisStn] = np.exp(qDaily[thisStn])

    nYears = int(qDaily.index.year.max() - qDaily.index.year.min())
    nSites = len(qDaily.columns)
    #print(qDaily.info())

    nR = min(1000, pow(nYears,2))
    num_realizations = [nYears, nR] #lall's criterion
    num_years = [nYears, 1]

    for r,y in zip(num_realizations,num_years) : #this is not about array index here
        print(r,'realisations -',y,'years')
        Qd_cg = base.combined_generator(qDaily, r, y, name=name)

        for i in Qd_cg.keys():
            for thisStn in Qd_cg[i].columns:
                if 'evap' in thisStn:
                    # back-transform evaporation
                    Qd_cg[i][thisStn] = np.log(Qd_cg[i][thisStn])
            Qd_cg[i].to_csv('../validation/synthetic/r'+str(i)+'X-y'+str(y)+'-daily.csv')
            monDf = convert_data_to_monthly(Qd_cg[i])
            for thisStn in monDf.columns:
                if 'evap' in thisStn:
                    # divide evaporation by 86400 (s/day) to get total monthly evap in mm/month
                    monDf[i][thisStn] /= 86400
            monDf.to_csv('../validation/synthetic/r'+str(i)+'X-y'+str(y)+'-monthly.csv')

if __name__ == "__main__":
    main()






