import numpy as np
import pandas as pd

qDaily = pd.csv_read('../data/Qdaily.txt')
nYears = qDaily.index.max() - qDaily.index.min()
nSites = len(qDaily.columns)
print(qDaily.info())

num_realizations = [100, 1000]
num_years = [100, 1]

for k in range(1:length(num_realizations)):
    Qd_cg = combined_generator(Qdaily, num_realizations[k], num_years[k] );
    """
    # back-transform evaporation
    Qd_cg[:,:,4] = log[Qd_cg[:,:,4]];

    # write simulated data to file
    for i=1:Nsites:
        q_ = [];
        for j=1:num_realizations(k)
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

def combined_generator(hist_data, nR, nY ) :
    # generation of monthly data via Kirsch et al. (2013):
    # Kirsch, B. R., G. W. Characklis, and H. B. Zeff (2013), 
    # Evaluating the impact of alternative hydro-climate scenarios on transfer 
    # agreements: Practical improvement for generating synthetic streamflows, 
    # Journal of Water Resources Planning and Management, 139(4), 396–406.
    QQg = monthly_main(hist_data, nR, nY );
    Qh = convert_data_to_monthly(hist_data);
    Nyears = size(Qh{1},1);

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

def monthly_main( hist_data, nR, nY ):
    # from daily to monthly
    #Qh = convert_data_to_monthly(hist_data) ;
    #use pandas transform here
    #create a monthly df with sums and adjusted from cumecs to cu.m over month
    Qh = df.groupby([lambda x: x.year, lambda x: x.month]).sum()
    Qh['Date'] = Qh.index.map(lambda x: date(x[0], x[1], 1))
    Qh.reset_index(drop=True, inplace=True)
    Qh.set_index('Date',inplace=True)
    Qh.index = pd.to_datetime(Qh.index)

    # generation
    for r in range(1:nR):
        Qs = monthly_gen(Qh, nY);
        for k=1:Nsites:
            qq{k}(r,:) = reshape(Qs{k}',1,[]);
        
    # output matrix
    Qgen = nan(nR,nY*12,Nsites) ;
    for k in range(1:Nsites):
        Qgen(:,:,k)= qq{k} ;

    return Qgen

def monthly_gen(Q_historical, num_years, p=None, n=None)
    nQ_historical = len(Q_historical);
        
    num_years = num_years+1; # this adjusts for the new corr technique
    if p == None and n == None:
        nQ = nQ_historical;
    else:
        n = n-1; # (input n=2 to double the frequency, i.e. repmat 1 additional time)
        nQ = nQ_historical + floor(p*nQ_historical+1)*n;
    
    Random_Matrix = randi(nQ, num_years, 12);
    
    for k in range(1:npoints):
        Q_matrix = Q_historical[k];
        """
        if  p != None and n != None:
            temp = sort(Q_matrix);
            app = temp(1:ceil(p*nQ_historical),:); # find lowest p# of values for each month
            Q_matrix = vertcat(Q_matrix, repmat(app, n, 1));
        """
        logQ = np.log10(Q_matrix);

        monthly_mean = [];
        monthly_stdev = [];
        Z = [];

        for i in range(1:13):
            monthly_mean.append(mean(logQ[i]))
            monthly_stdev.append(std(logQ[i]))
            Z.append((logQ[i] - monthly_mean[i]) / monthly_stdev[i])
        
        Z_vector = reshape(Z',1,[]);
        Z_shifted = reshape(Z_vector(7:(nQ*12-6)),12,[])';

        # The correlation matrices should use the historical Z's
        # (the "apped years" do not preserve correlation)
        U = chol_corr(Z(1:nQ_historical,:));
        U_shifted = chol_corr(Z_shifted(1:nQ_historical-1,:));

        for i=1:12
            Qs_uncorr(:,i) = Z(Random_Matrix(:,i), i);
        
        
        Qs_uncorr_vector = reshape(Qs_uncorr(:,:)',1,[]);
        Qs_uncorr_shifted(:,:) = reshape(Qs_uncorr_vector(7:(num_years*12-6)),12,[])';
        Qs_corr(:,:) = Qs_uncorr(:,:)*U;
        Qs_corr_shifted(:,:) = Qs_uncorr_shifted(:,:)*U_shifted;

        Qs_log(:,1:6) = Qs_corr_shifted(:,7:12);
        Qs_log(:,7:12) = Qs_corr(2:num_years,7:12);

        for i=1:12
            Qsk(:,i) = exp(Qs_log(:,i)*monthly_stdev(i) + monthly_mean(i));
        
        
        Qs{k} = Qsk;
    
    
    return Qs

