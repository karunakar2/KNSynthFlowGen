import glob

import numpy as np
import pandas as pd

from datetime import date

from base import *

fileList = glob.glob('../data/*.csv')
#print(fileList)

def writeExcelGrid(df,fileTag):
    #transform data and write to file in 'row wise years' and 'columnwise dates'
    for thisSite in df.columns:
        redDf = df[thisSite]
        temp = redDf.to_numpy()
        temp = temp.reshape(-1,365)
        np.savetxt('../validation/historical/'+thisSite+'-'+fileTag+'.csv', temp, delimiter=',', fmt='%1.4e')

df = pd.DataFrame()
for fName in fileList:
    thisDf = pd.read_csv(fName, parse_dates=["Date"])
    thisDf.set_index('Date',inplace=True)
    df = thisDf.copy() if len(df) < 1 else df.merge(thisDf, on='Date')
    #print(thisDf.head())

#the leap year extra day is not favourable and hence removed
df = df[~((df.index.month == 2) & (df.index.day == 29))]

#sort values by date
df.sort_values(by='Date', inplace=True)

#print(df.head())
#for further processing
df.to_csv('../data/Qdaily_.txt')

#for validation stuff
df.to_csv('../data/Qdaily.txt', header = False, index = False)
writeExcelGrid(df,'daily')

qMonthly = convert_data_to_monthly(df)
for thisStn in qMonthly.keys():
    if 'evap' in thisStn:
        qMonthly[thisStn] /= 24*3600 #cum/s to cum/day
    qMonthly[thisStn].to_csv('../validation/historical/'+thisStn+'-monthly.csv', header = False, index = False)


