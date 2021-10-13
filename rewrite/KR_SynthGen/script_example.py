import numpy as np
import pandas as pd

from base import *

##---------------------------------------------------------------------------------
def combined_generator(hist_data, nR, nY ) :
    # generation of monthly data via Kirsch et al. (2013):
    # Kirsch, B. R., G. W. Characklis, and H. B. Zeff (2013), 
    # Evaluating the impact of alternative hydro-climate scenarios on transfer 
    # agreements: Practical improvement for generating synthetic streamflows, 
    # Journal of Water Resources Planning and Management, 139(4), 396–406.
    QQg = monthly_main(hist_data, nR, nY );
    Qh = convert_data_to_monthly(hist_data);
    Nyears = len(Qh[1])

"""
    # disaggregation from monthly to daily time step as in Nowak et al. (2010):
    # Nowak, K., Prairie, J., Rajagopalan, B., & Lall, U. (2010). 
    # A nonparametric stochastic approach for multisite disaggregation of 
    # annual to daily streamflow. Water Resources Research, 46(8).

    # Find K-nearest neighbors (KNN) in terms of total monthly flow and 
    # randomly select one for disaggregation. Proportionally scale the flows in
    # the selected neighbor to match the synthetic monthly total. To
    # disaggregate Jan Flows, consider all historical January totals +/- 7
    # days, etc.
    Dt = 3600*24;
    DaysPerMonth = [31 28 31 30 31 30 31 31 30 31 30 31];
    D = zeros(nR,365*nY,Nsites);

    # concatenate last 7 days of last year before first 7 days of first year
    # and first 7 days of first year after last 7 days of last year
    nrows = size(hist_data,1);
    extra_hist_data = [hist_data(nrows-7:nrows,:); hist_data; hist_data(1:8,:)];

    # find monthly totals for all months +/- 7 days
    for i=1:12
        count = 1;
        if i == 1 || i == 12
            nTotals = Nyears*15-7; # 7 less shifts in first and last month
        else
            nTotals = Nyears*15;
        Qmonthly_shifted = zeros(nTotals,Nsites);
        indices = zeros(nTotals,2);
        for k = 1:15
            shifted_hist_data = extra_hist_data(k:k+nrows-1,:);
            Qh = convert_data_to_monthly(shifted_hist_data);
            for j=1:Nsites
                if i == 1 && k<8
                    Qh{j} = Qh{j}(2:size(Qh{j},1),i); # remove first year
                elseif i == 12 && k>8
                    Qh{j} = Qh{j}(1:(size(Qh{j},1)-1),i); # remove last year
                Qmonthly_shifted(count:(count+size(Qh{j},1)-1),j) = Qh{j}(:,1);
            if i == 1 && k<8
                indices(count:(count+size(Qh{j},1)-1),1) = 2:(size(Qh{j},1)+1);
            else
                indices(count:(count+size(Qh{j},1)-1),1) = 1:size(Qh{j},1);
            indices(count:(count+size(Qh{j},1)-1),2) = k;
            count = count + size(Qh{j},1);
        Qtotals{i} = Qmonthly_shifted;
        Qindices{i} = indices;

    for r=1:nR
        dd = [];
        for i=1:nY*12
            #monthly value for all sites
            Z = QQg(r,i,:);
            #KNN and weights
            month = mod(i,12);
            if month == 0
                month = 12;
            [KNN_id, W] = KNN_identification(Z, Qtotals, month);
            Wcum = cumsum(W);
                
            #sampling of one KNN
            py = KNN_sampling(KNN_id, Qindices{month}, Wcum, hist_data, month);
            d = zeros(Nsites,DaysPerMonth(month));
            for j=1:Nsites
                d(j,:) = py(:,j).*Z(1,1,j);
            dd = [dd d];
        
        D(r,:,:) = dd'/Dt;
"""
##---------------------------------------------------------------------------------


qDaily = pd.read_csv('../data/Qdaily.txt', parse_dates=["Date"])
qDaily.set_index('Date',inplace=True)
for thisStn in qDaily.columns:
    if 'evap' in thisStn:
        qDaily[thisStn] = np.exp(qDaily[thisStn])

nYears = int(qDaily.index.year.max() - qDaily.index.year.min())
nSites = len(qDaily.columns)
#print(qDaily.info())

num_realizations = [pow(nYears,2), 1000] #lall's criterion
num_years = [pow(nYears,2), 1]

for k in range(1,len(num_realizations)):
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








