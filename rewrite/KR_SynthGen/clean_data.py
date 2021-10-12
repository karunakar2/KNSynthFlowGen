import glob

import numpy as np
import pandas as pd

from datetime import date

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

#create a monthly df with sums and adjusted from cumecs to cu.m over month
mDf = df.groupby([lambda x: x.year, lambda x: x.month]).sum()
mDf['Date'] = mDf.index.map(lambda x: date(x[0], x[1], 1))
mDf.reset_index(drop=True, inplace=True)
mDf.set_index('Date',inplace=True)
mDf.index = pd.to_datetime(mDf.index)
#print(mDf.info())

#conversion of flow to cum/day
for thisCol in mDf.columns:
    if 'evap' not in thisCol:
        mDf[thisCol] *= 24*3600 #cum/s to cum/day

writeExcelGrid(mDf,'monthly')
