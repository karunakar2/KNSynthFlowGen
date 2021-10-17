import numpy as np
import pandas as pd

from datetime import date

#Validated functions
#Sum aggregates the daily data to monthly values
def convert_data_to_monthly(df:pd.DataFrame()) ->{'station name':pd.DataFrame()}:
    mDf = df.groupby([lambda x: x.year, lambda x: x.month]).sum()
    #print(mDf.head())
    """
    mDf = df.groupby([lambda x: x.year, lambda x: x.month]).sum()
    mDf['Date'] = mDf.index.map(lambda x: date(x[0], x[1], 1))
    mDf.reset_index(drop=True, inplace=True)
    mDf.set_index('Date',inplace=True)
    mDf.index = pd.to_datetime(mDf.index)
    """

    yearList = list(mDf.index.get_level_values(0))
    monthList = list(range(1,13))
    qMonthly = {}
    #inititalisation in pandas, defeats the purpose of being dynamic
    for thisSite in mDf.columns:
        qMonthly[thisSite] = pd.DataFrame() #columns=monthList) #, index=yearList)
        #populate the dataframes
        for thisYear in range(min(yearList), max(yearList)):
            redDf = mDf.loc[thisYear,thisSite]
            data = {(key+1):value for key,value in enumerate(redDf.values)}
            tempDf = pd.DataFrame(data= data, columns=monthList, index=[thisYear])
            #print(tempDf.head())
            qMonthly[thisSite] = qMonthly[thisSite].append(tempDf)
        #print(qMonthly[thisSite].head())
    return qMonthly

#patches the data in front and end (N+D+N days) to sift through adjustments
def patchData(hist_data:pd.DataFrame(), period:int, debug:bool=0):
    preDf = hist_data.tail(period).copy().reset_index(drop=True)
    preDf['Date'] = pd.date_range(hist_data.index.min(), freq='-1D', periods=period)
    preDf.set_index('Date',inplace=True)
    preDf.sort_values(by=['Date'],inplace=True)
    posDf = hist_data.head(period).copy().reset_index(drop=True)
    posDf['Date'] = pd.date_range(hist_data.index.max(), periods=period)
    posDf.set_index('Date',inplace=True)
    extra_hist_data = pd.concat([preDf,hist_data,posDf])
    if debug != 0:
        print(extra_hist_data.head())
        print(extra_hist_data.tail())
        print(len(hist_data),len(extra_hist_data))
    
    return extra_hist_data

##-------------------------------------------------------------------------------
##Unverified functions
eps = np.finfo(float).eps
def chol_corr(Z):
# compute cholesky decomp of correlation matrix of columns of Z
# attempts to repair non-positive-definite matrices
# http://www.mathworks.com/matlabcentral/answers/6057-repair-non-positive-definite-correlation-matrix
# rank-1 update followed by rescaling to get unit diagonal entries
    U = None
    R = np.corrcoef(Z)
    try :
        U = np.linalg.cholesky(R)
    #will raise LinAlgError if the matrix is not positive definite.
    #https://stackoverflow.com/a/16266736
    except Exception as er:
        if er == 'LinAlgError':
            posDef = False
            while posDef:
                try :
                    k = min([min(np.real(np.linalg.eig(R))) -1*eps])
                    R = R - k*np.eye(R.size)
                    R = R/R[0,0] #adjusted for indices
                    U = np.linalg.cholesky(R)
                    posDef = False
                except Exception as err:
                    if err == 'LinAlgError':
                        posDef = True
                    else:
                        print(err)
            else:
                print(er)
    return U

def monthly_main( hist_data:pd.DataFrame(), nR:int, nY:int ):
    # from daily to monthly
    Qh_mon = convert_data_to_monthly(hist_data) 
    Nsites = len(Qh)

    ##qq ={}
    #qq = np.empty((nR,nY*12,Nsites)).fill(np.nan)
    """
    for r in range(0,nR):
        Qs = monthly_gen(Qh, nY)
        for k in range(1,Nsites):
            #print(len(Qs[k]),len(Qs[k][1]))
            #qq{k}(r,:) = reshape(Qs{k}',1,[])
            temp = np.ravel(Qs[k].T) #Qs[k] should be 12 elements but returns 816 i.e. (12months *68Years) here
            try:
                len(qq[k])
            except KeyError:
                qq[k] = np.empty((nR,len(temp)))
            qq[k][r,:] = temp
    
    # output matrix
    #Qgen = nan(nR,nY*12,Nsites) 
    Qgen = np.empty((nR,nY*12,Nsites)) #.fill(np.nan)
    for k in range(1,Nsites):
        Qgen[:,:,k]= qq[k] 
    return Qgen
    """
    # output matrix
    #Qgen = nan(nR,nY*12,Nsites) 
    ##Qgen = np.empty((nR,nY*12,Nsites)) #.fill(np.nan)
    Qgen = {}
    for r in range(0,nR):
        Qmon_gen = monthly_gen(Qh_mon, nY) #formerly Qs
        """
        temp = []
        for k in range(0,Nsites):
            #print(len(Qs[k]),len(Qs[k][1]))
            #qq{k}(r,:) = reshape(Qs{k}',1,[])
            temp.append(np.ravel(Qs[k].T))
        temp = np.array(temp).T
        #print(temp.shape)
        Qgen[r] = temp
        #Qgen[r,:,k] = temp.reshape((1,temp.size,1))
        """
        #why do above transformations
        Qgen[k] = Qmon_gen
    
    #print(Qgen,'Qgen')
    ##print(len(Qgen),Qgen[1].shape)
    return Qgen

def monthly_gen(q_historical, num_years, p=None, n=None):
    """
    q_historical is monthly
    """
    qH_stns = list(q_historical.keys()) #keys are station names
    nPoints = len(qH_stns) #no of donor stations
    nQ_historical = len(q_historical[qH_stns[0]]) #no of years(with months) of data
    #future - check all dfs are of same size in qhistorical
    
    num_years = num_years+1 # this adjusts for the new corr technique
    if p == None and n == None:
        nQ = nQ_historical
    else:
        n = n-1 # (input n=2 to double the frequency, i.e. repmat 1 additional time)
        nQ = nQ_historical + math.floor(p*nQ_historical+1)*n
    
    #Random_Matrix = randi(nQ, num_years, 12)
    #Random_Matrix = np.random.randint(num_years, size=(nQ,12))
    Random_Matrix = np.random.randint(nQ, size=(num_years,12))
    #print(Random_Matrix.shape,'rmsize')
    
    Qs = {}
    for k in q_historical.keys():
        Q_matrix = q_historical[k]
        #print(Q_matrix.shape)
        """
        if  p != None and n != None:
            temp = sort(Q_matrix)
            app = temp(1:ceil(p*nQ_historical),:) # find lowest p# of values for each month
            Q_matrix = vertcat(Q_matrix, repmat(app, n, 1))
        """
        logQ = Q_matrix.apply(lambda x: np.log(x))
        monthly_mean = []
        monthly_stdev = []
        Z = []
        for mon in logQ.columns: 
            m_series = logQ[mon].values
            m_mean = np.mean(m_series)
            m_stdev = np.std(m_series)
            Z.append((m_series - m_mean) / m_stdev) #zscore here and double dim due to logQ[i] series
            monthly_mean.append(m_mean)
            monthly_stdev.append(m_stdev)
        Z = np.array(Z) #list of lists to array format
        #print(Z.shape)
        
        #Z_vector = reshape(Z',1,[])
        Z_vector = np.ravel(Z.T)
        Z_JulyToJune = Z_vector[6:-6] #zero index so starts at 6, uptill last 6 months
        #Z_shifted = reshape(Z_vector(7:(nQ*12-6)),12,[])' #july of start year to june of end year
        Z_shifted = np.transpose(Z_JulyToJune.reshape(-1,12))
        #print(Z_shifted.shape,'zshiftshape')
        
        # The correlation matrices should use the historical Z's
        # (the "appendedd years" do not preserve correlation)
        U = chol_corr(Z[:nQ_historical-1,:]) #index adjusted for python
        U_shifted = chol_corr(Z_shifted[:nQ_historical-2,:])
        #print('b',U.shape)
        
        Qs_uncorr = []
        for i in range(0,12): #python index related
            #print('z RM',Z.shape,Random_Matrix.shape)
            #Qs_uncorr[:,i] = Z[Random_Matrix[:,i], i]
            ##Qs_uncorr.append(Z[Random_Matrix[:,i],i])
            Qs_uncorr.append(Z[i,Random_Matrix[:,i]])
        Qs_uncorr = np.array(Qs_uncorr).T
        #Qs_corr(:,:) = Qs_uncorr(:,:)*U
        #print('c',Qs_uncorr.shape)
        Qs_corr = np.dot(Qs_uncorr,U)
        
        #Qs_uncorr_vector = reshape(Qs_uncorr(:,:)',1,[])
        Qs_uncorr_vector = np.ravel(Qs_uncorr)
        #Qs_uncorr_shifted(:,:) = reshape(Qs_uncorr_vector(7:(num_years*12-6)),12,[])'
        Qs_uncorr_vector_JJ = Qs_uncorr_vector[6:-6]
        Qs_uncorr_shifted = Qs_uncorr_vector_JJ.reshape(-1,12)
        #Qs_corr_shifted(:,:) = Qs_uncorr_shifted(:,:)*U_shifted
        Qs_corr_shifted = np.dot(Qs_uncorr_shifted,U_shifted)

        Qs_log = np.empty((len( Qs_corr_shifted ),len( Qs_corr_shifted [0])))
        #print(Qs_log)
        Qs_log[:,0:5] = Qs_corr_shifted[:,6:11] #indices adjusted
        #Qs_log[:,6:11] = Qs_corr[1:num_years-1,6:11]
        Qs_log[:,6:11] = Qs_corr[:-1,6:11] #not sure should we exclude last or first?
        #print(Qs_log.shape)

        Qsk = np.empty((len( Qs_log ),len( Qs_log[0])))
        for i in range(0,12):
            Qsk[:,i] = np.exp(Qs_log[:,i]*monthly_stdev[i] + monthly_mean[i])
        QskDf = pd.DataFrame(Qsk, columns = Q_matrix.columns, index = Q_matrix.index)
        Qs[k] = QskDf
    
    #print(Qs)
    return Qs

###---------------------yet to finish this part -----------------------------------

def KNN_identification( Z, Qtotals, month, K=None ):
    # [KNN_id, W] = KNN_identification( Z, Qtotals, month, k )
    #
    # Identification of K-nearest neighbors of Z in the historical annual data
    # z and computation of the associated weights W.
    #
    # Input:    Z = synthetic datum (scalar)
    #           Qtotals = total monthly flows at all sites for all historical months 
    #             within +/- 7 days of the month being disaggregated
    #           month = month being disaggregated
    #           k = number of nearest neighbors (by default k=n_year^0.5
    #             according to Lall and Sharma (1996))
    # Output:   KNN_id = indices of the first K-nearest neighbors of Z in the
    #             the historical annual data z
    #           W = nearest neighbors weights, according to Lall and Sharma
    #             (1996): W(i) = (1/i) / (sum(1/i)) 
    #
    # MatteoG 31/05/2013

    # Ntotals is the number of historical monthly patterns used for disaggregation.
    # A pattern is a sequence of ndays of daily flows, where ndays is the
    # number of days in the month being disaggregated. Patterns are all
    # historical sequences of length ndays beginning within 7 days before or
    # after the 1st day of the month being disaggregated.
    Ntotals = len(Qtotals[month])
    if K == None:
        K = round(sqrt(Ntotals))

    # nearest neighbors identification
    # only look at neighbors from the same month +/- 7 days
    Nsites = len(Qtotals[month][0])
    delta = np.empty((Ntotals,1)) # first and last month have 7 less possible shifts
    for i in range(1,Ntotals):
        for j in range(1,Nsites):
                delta[i] += (Qtotals[month][i,j]-Z[1,1,j])^2

    #Y = [[1:size(delta,1)]', delta ] 
    Y = pd.series(np.array(delta))
    #Y_ord = sortrows(Y, 2)
    Y_ord = Y.sort_values()
    ###What to do here??
    KNN_id = Y_ord[1:K,1] 

    # computation of the weights
    f = np.array(range(1,K))
    f1 = 1/f
    W = f1 / sum(f1) 

    return KNN_id, W

def KNN_sampling( KNN_id, indices, Wcum, Qdaily, month ):
    # py = KNN_sampling( KKN_id, indices, Wcum, Qdaily, month )
    #
    # Selection of one KNN according to the probability distribution defined by
    # the weights W.
    #
    # Input:    KNN_id = indices of the first K-nearest neighbors
    #           indices = n x 2 matrix where n is the number of monthly totals
    #             and the 2 columns store the historical year in which each
    #             monthly total begins, and the number of shift index
    #             where 1 is 7 days earlier and 15 is 7 days later
    #           Wcum = cumulated probability for each nearest neighbor
    #           Qdaily = historical data
    #           month = month being disaggregated
    # Output:   py = selected proportion vector corresponding to the sampled
    #             shifted historical month
    #           yearID = randomly selected monthly total (row to select from indices)
    #
    # MatteoG 31/05/2013

    #Randomly select one of the k-NN using the Lall and Sharma density estimator
    #r = rand 
    r = np.random.rand(1,1)
    Wcum = [0, Wcum]  #not sure whats happening here
    for i in range(0,len(Wcum)) : #python index adjusted
        if (r > Wcum[i]) and (r <= Wcum[i+1]) :
            KNNs = i 
    yearID = KNN_id[KNNs]

    # concatenate last 7 days of last year before first 7 days of first year
    # and first 7 days of first year after last 7 days of last year
    nrows = len(Qdaily)
    ##Qdaily = [Qdaily[nrows-7:nrows,:] Qdaily Qdaily[1:8,:]]
    Qdaily = patchData(Qdaily,7)
    
    #well this needs fixing
    """
    # shift historical data to get nearest neighbor corresponding to yearID
    year = indices(yearID,1) #find the year that best fits
    k = indices(yearID,2) #+/- 7 day shift range that has high Weightage
    shifted_Qdaily = Qdaily(k:k+nrows-1,:)

    DaysPerMonth = [31 28 31 30 31 30 31 31 30 31 30 31]
    start = 365*(year-1) + sum(DaysPerMonth(1:(month-1)))+1
    dailyFlows = shifted_Qdaily(start:start+DaysPerMonth(month)-1,:)

    py = zeros(size(dailyFlows))
    for i=1:size(Qdaily,2)
        py(:,i) = dailyFlows(:,i)/sum(dailyFlows(:,i))

    return [py, yearID] 
    """
    #year = ##
    #k = #offset
    shifted_Qdaily = Qdaily.iloc[k:k+nrows,:]
    dailyFlows = shifted_Qdaily[shifted_Qdaily.index.year== year and shifted_Qdaily.index.month == month].copy()
    
    for thisColumn in dailyFlows.columns:
        py[:,thisColumn] = dailyFlows[thisColumn]/sum(dailyFlows[thisColumn])
        
    return [py, yearID]


##---------------------------------------------------------------------------------
DaysPerMonth = [31,28,31,30,31,30,31,31,30,31,30,31]
def combined_generator(hist_data, nR, nY ) :
    
    # generation of monthly data via Kirsch et al. (2013):
    # Kirsch, B. R., G. W. Characklis, and H. B. Zeff (2013), 
    # Evaluating the impact of alternative hydro-climate scenarios on transfer 
    # agreements: Practical improvement for generating synthetic streamflows, 
    # Journal of Water Resources Planning and Management, 139(4), 396–406.
    Q_Qg = monthly_main(hist_data, nR, nY )
    Qh = convert_data_to_monthly(hist_data)
    Qh_stns = list(Qh.keys())
    Nyears = len(Qh[Qh_stns[0]])

    # disaggregation from monthly to daily time step as in Nowak et al. (2010):
    # Nowak, K., Prairie, J., Rajagopalan, B., & Lall, U. (2010). 
    # A nonparametric stochastic approach for multisite disaggregation of 
    # annual to daily streamflow. Water Resources Research, 46(8).

    # Find K-nearest neighbors (KNN) in terms of total monthly flow and 
    # randomly select one for disaggregation. Proportionally scale the flows in
    # the selected neighbor to match the synthetic monthly total. To
    # disaggregate Jan Flows, consider all historical January totals +/- 7
    # days, etc.
    ###Dt = 3600*24
    ###DaysPerMonth = [31 28 31 30 31 30 31 31 30 31 30 31]
    ###D = np.empty((nR,365*nY,Nsites))

    # concatenate last 7 days of last year before first 7 days of first year
    # and first 7 days of first year after last 7 days of last year
    
    #extra_hist_data = [hist_data[nrows-7:nrows,:] hist_data hist_data(1:8,:)]
    extra_hist_data = patchData(hist_data,7)
    
"""
    # find monthly totals for all months +/- 7 days
    for i=1:12
        count = 1
        if i == 1 || i == 12
            nTotals = Nyears*15-7 # 7 less shifts in first and last month
        else
            nTotals = Nyears*15
        Qmonthly_shifted = zeros(nTotals,Nsites)
        indices = zeros(nTotals,2)
        for k = 1:15
            shifted_hist_data = extra_hist_data(k:k+nrows-1,:)
            Qh = convert_data_to_monthly(shifted_hist_data)
            for j=1:Nsites
                if i == 1 && k<8
                    Qh{j} = Qh{j}(2:size(Qh{j},1),i) # remove first year
                elseif i == 12 && k>8
                    Qh{j} = Qh{j}(1:(size(Qh{j},1)-1),i) # remove last year
                Qmonthly_shifted(count:(count+size(Qh{j},1)-1),j) = Qh{j}(:,1)
            if i == 1 && k<8
                indices(count:(count+size(Qh{j},1)-1),1) = 2:(size(Qh{j},1)+1)
            else
                indices(count:(count+size(Qh{j},1)-1),1) = 1:size(Qh{j},1)
            indices(count:(count+size(Qh{j},1)-1),2) = k
            count = count + size(Qh{j},1)
        Qtotals{i} = Qmonthly_shifted
        Qindices{i} = indices
    """
    nrows = len(hist_data)
    Qtotals = {}
    for k in range(1,15):
        Qmonthly_shifted = (Qh.iloc[k:k+nrows,:]).copy()
        Qmonthly_shifted.reset_index(inplace=True,drop=True)
        #Qmonthly_shifted.index = Qh.index
        Qmonthly_shifted.reindex_like(Qh)
        Qtotals[k] = Qmonthly_shifted
        
    #what should we do about indices here?
    
    """
    for r=1:nR
        dd = []
        for i=1:nY*12
            #monthly value for all sites
            Z = Q_Qg[r,i,:]
            #KNN and weights
            month = mod(i,12)
            if month == 0
                month = 12
            [KNN_id, W] = KNN_identification(Z, Qtotals, month)
            Wcum = cumsum(W)
                
            #sampling of one KNN
            py = KNN_sampling(KNN_id, Qindices{month}, Wcum, hist_data, month)
            d = zeros(Nsites,DaysPerMonth(month))
            for j=1:Nsites
                d(j,:) = py(:,j).*Z(1,1,j)
            dd = [dd d]
        
        D(r,:,:) = dd'/Dt
"""
    for r in range(0,nR):
    dd = []
    for i, month in zip(range(0,nY*12), range(0,12)):
        z = Q_Qg[r,i,:]
        [KNN_id, W] = KNN_identification(Z, Qtotals, month)
        py, _ = KNN_sampling(KNN_id, Qindices{month}, Wcum, hist_data, month)
        #d = zeros(Nsites,DaysPerMonth(month))
        d = np.empty((Nsites,DaysPerMonth(month)))
        for j in range(0:Nsites):
            d[j,:] = py[:,j].*Z[1,1,j]
        dd = [dd d] #what is this here? append?
        
        D[r,:,:] = dd'/Dt #is this compliment or ?
    
    return D
##---------------------------------------------------------------------------------
