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
        #print(Qs,'Qs')
        for k in range(1,Nsites):
            #qq{k}(r,:) = reshape(Qs{k}',1,[]);
            qq[k][r,:] = np.ravel(Qs[k].T)
    
    # output matrix
    #Qgen = nan(nR,nY*12,Nsites) ;
    Qgen = np.empty((nR,nY*12,Nsites)).fill(np.nan)
    
    for k in range(0,Nsites):
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


