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
        filePtr = open('../validation/historical/'+thisSite+'-'+fileTag+'.csv','w')
        for thisYear in range(df.index.year.min(),df.index.year.max()+1):
            filePtr.write((',').join(df[df.index.year == thisYear][thisSite].values.astype(str)))
            filePtr.write('\n')
    filePtr.close()

df = pd.DataFrame()
for fName in fileList:
    thisDf = pd.read_csv(fName, parse_dates=["Date"])
    thisDf.set_index('Date',inplace=True)
    if len(df) < 1:
        df = thisDf.copy()
    else:
        df = df.merge(thisDf, on='Date')
    #print(thisDf.head())

#the leap year extra day is not favourable and hence removed
df = df[~((df.index.month == 2) & (df.index.day == 29))]

#print(df.head())
df.to_csv('../data/Qdaily.txt')
writeExcelGrid(df,'daily')

qMonthly = convert_data_to_monthly(df)
for thisStn in qMonthly.keys():
    if 'evap' not in thisStn:
        qMonthly[thisStn] *= 24*3600 #cum/s to cum/day
    qMonthly[thisStn].to_csv('../validation/historical/'+thisStn+'-monthly.csv')


