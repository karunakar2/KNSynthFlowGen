import numpy as np
import pandas as pd

from base import *


qDaily = pd.read_csv('../data/Qdaily.txt', parse_dates=["Date"])
qDaily.set_index('Date',inplace=True)
for thisStn in qDaily.columns:
    if 'evap' in thisStn:
        qDaily[thisStn] = np.exp(qDaily[thisStn])

nYears = int(qDaily.index.year.max() - qDaily.index.year.min())
nSites = len(qDaily.columns)
#print(qDaily.info())

nR = min(1000, pow(nYears,2))
num_realizations = [nR, nYears] #lall's criterion
num_years = [nR, 1]

for k in range(0,len(num_realizations)):
    Qd_cg = combined_generator(qDaily, num_realizations[k], num_years[k] )
    """
    # back-transform evaporation
    Qd_cg[:,:,4] = log[Qd_cg[:,:,4]];

    # write simulated data to file
    for i in range(1:Nsites):
        q_ = [];
        for j in range(1:num_realizations(k)):
            qi = nan(365*num_years(k),1);
            qi(1:size(Qd_cg,2)) = Qd_cg(j,:,i)';
            q_ = [q_ reshape(qi,365,num_years(k))];
        
        Qd2(:,i) = reshape(q_(:),[],1);
        saveQ = reshape(Qd2(:,i), num_years(k)*365, num_realizations(k))';
        dlmwrite(['./../validation/synthetic/' sites{i} dimensions{k} '-daily.csv'], saveQ);
    synMonthlyQ = convert_data_to_monthly(Qd2);
    # divide evaporation by 86400 (s/day) to get total monthly evap in mm/month
    synMonthlyQ{4} = synMonthlyQ{4}/86400;
    for i=1:Nsites
        saveMonthlyQ = reshape(synMonthlyQ{i}',12*num_years(k),num_realizations(k))';
        dlmwrite(['./../validation/synthetic/' sites{i} dimensions{k} '-monthly.csv'], saveMonthlyQ);
    dlmwrite(['./../validation/synthetic/Qdaily' dimensions{k} '.csv'], Qd2);
    clear Qd2;
    """








