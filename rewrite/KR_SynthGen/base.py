import numpy as np
import pandas as pd

from datetime import date

#Validated functions
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

##-------------------------------------------------------------------------------
##Unverified functions
eps = np.finfo(float).eps
def chol_corr(Z):
# compute cholesky decomp of correlation matrix of columns of Z
# attempts to repair non-positive-definite matrices
# http://www.mathworks.com/matlabcentral/answers/6057-repair-non-positive-definite-correlation-matrix
# rank-1 update followed by rescaling to get unit diagonal entries
    U = None
    R = np.corrcoef(Z);
    try :
        U = np.linalg.cholesky(R);
    #will raise LinAlgError if the matrix is not positive definite.
    #https://stackoverflow.com/a/16266736
    except Exception as er:
        if er == 'LinAlgError':
            posDef = False
            while posDef:
                try :
                    k = min([min(np.real(np.linalg.eig(R))) -1*eps]);
                    R = R - k*np.eye(R.size);
                    R = R/R[1,1];
                    U = np.linalg.cholesky(R);
                    posDef = False
                except Exception as err:
                    if err == 'LinAlgError':
                        posDef = True
                    else:
                        print(err)
            else:
                print(er)
    return U

def monthly_main( hist_data, nR, nY ):
    # from daily to monthly
    Qh = convert_data_to_monthly(hist_data) ;
    Nsites = len(Qh)

    qq ={}
    #qq = np.empty((nR,nY*12,Nsites)).fill(np.nan)
    # generation
    for r in range(0,nR):
        Qs = monthly_gen(Qh, nY);
        for k in range(1,Nsites):
            #print(Qs[k],'Qs')
            print(len(Qs[k]),len(Qs[k][1]))
            raise Exception('stop here for me')
            #qq{k}(r,:) = reshape(Qs{k}',1,[]);
            temp = np.ravel(Qs[k].T) #Qs[k] should be 12 elements but returns 816 i.e. (12months *68Years) here
            try:
                len(qq[k])
            except KeyError:
                qq[k] = np.empty((nR,len(temp)))
            qq[k][r,:] = temp
    
    # output matrix
    #Qgen = nan(nR,nY*12,Nsites) ;
    Qgen = np.empty((nR,nY*12,Nsites)) #.fill(np.nan)
    for k in range(1,Nsites):
        Qgen[:,:,k]= qq[k] ;
    return Qgen

def monthly_gen(q_historical, num_years, p=None, n=None):
    q_h_keys = list(q_historical.keys())
    nPoints = len(q_h_keys ) #no of donor stations
    nQ_historical = len(q_historical[q_h_keys[1]]); #no of years of data
    #future - check all dfs are of same size in qhistorical
    
    num_years = num_years+1; # this adjusts for the new corr technique
    if p == None and n == None:
        nQ = nQ_historical;
    else:
        n = n-1; # (input n=2 to double the frequency, i.e. repmat 1 additional time)
        nQ = nQ_historical + math.floor(p*nQ_historical+1)*n;
    
    #Random_Matrix = randi(nQ, num_years, 12);
    Random_Matrix = np.random.randint(num_years, size=(nQ,12))
    #print(Random_Matrix)
    
    Qs = {}
    for k in range(1,nPoints):
        Q_matrix = q_historical[q_h_keys[k]]
        #print(Q_matrix)
        """
        if  p != None and n != None:
            temp = sort(Q_matrix);
            app = temp(1:ceil(p*nQ_historical),:); # find lowest p# of values for each month
            Q_matrix = vertcat(Q_matrix, repmat(app, n, 1));
        """
        logQ = Q_matrix.apply(lambda x: np.log(x))
        monthly_mean = [];
        monthly_stdev = [];
        Z = [];
        for i in range(1,13):
            m_series = logQ[i].values
            m_mean = np.mean(m_series)
            m_stdev = np.std(m_series)
            Z.append((m_series - m_mean) / m_stdev) #zscore here and double dim due to logQ[i] series
            monthly_mean.append(m_mean)
            monthly_stdev.append(m_stdev)
        Z = np.array(Z) #list of lists to array format
        
        #Z_vector = reshape(Z',1,[]);
        Z_vector = np.ravel(np.array(Z))
        Z_JulyToJune = Z_vector[6:-6] #zero index so starts at 6, uptill last 6 months
        #Z_shifted = reshape(Z_vector(7:(nQ*12-6)),12,[])'; #july of start year to june of end year
        Z_shifted = np.transpose(Z_JulyToJune.reshape(-1,12))

        # The correlation matrices should use the historical Z's
        # (the "appendedd years" do not preserve correlation)
        U = chol_corr(Z[:nQ_historical-1,:]) #index adjusted for python
        U_shifted = chol_corr(Z_shifted[:nQ_historical-2,:])

        Qs_uncorr = []
        for i in range(0,12): #python index related
            #Qs_uncorr[:,i] = Z[Random_Matrix[:,i], i];
            Qs_uncorr.append(Z[Random_Matrix[:,i],i])
        Qs_uncorr = np.array(Qs_uncorr).T
        #Qs_corr(:,:) = Qs_uncorr(:,:)*U;
        Qs_corr = np.dot(Qs_uncorr,U)
        
        #Qs_uncorr_vector = reshape(Qs_uncorr(:,:)',1,[]);
        Qs_uncorr_vector = np.ravel(Qs_uncorr)
        #Qs_uncorr_shifted(:,:) = reshape(Qs_uncorr_vector(7:(num_years*12-6)),12,[])';
        Qs_uncorr_vector_JJ = Qs_uncorr_vector[6:-6]
        Qs_uncorr_shifted = Qs_uncorr_vector_JJ.reshape(-1,12)
        #Qs_corr_shifted(:,:) = Qs_uncorr_shifted(:,:)*U_shifted;
        Qs_corr_shifted = np.dot(Qs_uncorr_shifted,U_shifted)

        Qs_log = np.empty((len( Qs_corr_shifted ),len( Qs_corr_shifted [1])))
        #print(Qs_log)
        Qs_log[:,0:5] = Qs_corr_shifted[:,6:11] #indices adjusted
        #Qs_log[:,6:11] = Qs_corr[1:num_years-1,6:11]
        Qs_log[:,6:11] = Qs_corr[:-1,6:11] #not sure should we exclude last or first?

        Qsk = np.empty((len( Qs_log ),len( Qs_log[1])))
        for i in range(0,12):
            Qsk[:,i] = np.exp(Qs_log[:,i]*monthly_stdev[i] + monthly_mean[i]);
        
        Qs[k] = Qsk;
        
    return Qs

###---------------------yet to finish this part -----------------------------------

def KNN_identification( Z, Qtotals, month, K=None )
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
    Ntotals = len(Qtotals[month]);
    if K == None:
        K = round(sqrt(Ntotals))

    # nearest neighbors identification;
    # only look at neighbors from the same month +/- 7 days
    Nsites = len(Qtotals[month][1]);
    delta = np.empty((Ntotals,1)); # first and last month have 7 less possible shifts
    for i in range(1:Ntotals):
        for j in range(1:Nsites):
                delta[i] += (Qtotals[month][i,j]-Z[1,1,j])^2;

    #Y = [[1:size(delta,1)]', delta ] ;
    Y = [np.array(list(range(1,len(delta)))).T,np.array(delta)]
    #Y_ord = sortrows(Y, 2);
    ###What to do here??
    KNN_id = Y_ord[1:K,1] ;

    # computation of the weights
    f = np.array(range(1,K))
    f1 = 1./f;
    W = f1 ./ sum(f1) ;

    return KNN_id, W

def KNN_sampling( KNN_id, indices, Wcum, Qdaily, month )
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
    r = rand ;
    Wcum = [0, Wcum] ;
    for i = 1:len(Wcum)-1 :
        if (r > Wcum[i]) && (r <= Wcum[i+1]) :
            KNNs = i ;
    yearID = KNN_id( KNNs ) ;

    # concatenate last 7 days of last year before first 7 days of first year
    # and first 7 days of first year after last 7 days of last year
    nrows = len(Qdaily)
    Qdaily = [Qdaily[nrows-7:nrows,:]; Qdaily; Qdaily[1:8,:]];
    #what is this here, multitude of columns?
"""
    # shift historical data to get nearest neighbor corresponding to yearID
    year = indices(yearID,1);
    k = indices(yearID,2);
    shifted_Qdaily = Qdaily(k:k+nrows-1,:);

    DaysPerMonth = [31 28 31 30 31 30 31 31 30 31 30 31];
    start = 365*(year-1) + sum(DaysPerMonth(1:(month-1)))+1;
    dailyFlows = shifted_Qdaily(start:start+DaysPerMonth(month)-1,:);

    py = zeros(size(dailyFlows));
    for i=1:size(Qdaily,2)
        py(:,i) = dailyFlows(:,i)/sum(dailyFlows(:,i));

    return [py, yearID] 
"""
