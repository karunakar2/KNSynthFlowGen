#import sys
import warnings

import math
import numpy as np
import pandas as pd

from datetime import date

import shelve

import logging
logging.basicConfig(filename='./test.log', level=logging.DEBUG, format='%(asctime)s | %(levelname)s | %(funcName)s |%(message)s')


#Validated functions
#Sum aggregates the daily data to monthly values
def convert_data_to_monthly(df:pd.DataFrame(), debug=False) ->{'station name':pd.DataFrame()}:
    mDf = df.groupby([lambda x: x.year, lambda x: x.month]).sum() * 86400 
    #convert flow/s to flow/day * monthsum
    if debug:
        logging.info('mDf',mDf.head())
    """
    mDf = df.groupby([lambda x: x.year, lambda x: x.month]).sum()
    mDf['Date'] = mDf.index.map(lambda x: date(x[0], x[1], 1))
    mDf.reset_index(drop=True, inplace=True)
    mDf.set_index('Date',inplace=True)
    mDf.index = pd.to_datetime(mDf.index)
    """

    yearList = list(mDf.index.get_level_values(0).unique())
    if debug:
        logging.info('yearList',yearList)
    monthList = list(range(1,13))
    qMonthly = {}
    #inititalisation in pandas, defeats the purpose of being dynamic
    for thisSite in mDf.columns:
        qMonthly[thisSite] = pd.DataFrame()
        #populate the dataframes
        for thisYear in yearList: #range(min(yearList), max(yearList)+1):
            redDf = mDf.loc[thisYear,thisSite]
            data = {(key+1):value for key,value in enumerate(redDf.values)} 
            #key start at 0, align it to 1 i.e jan
            tempDf = pd.DataFrame(data= data, columns=monthList, index=[thisYear])
            if debug:
                logging.info('tempDf',tempDf.head())
            qMonthly[thisSite] = qMonthly[thisSite].append(tempDf)
        if debug:
            logging.info('qMonthly',thisSite,qMonthly[thisSite].head())
    return qMonthly

#patches the data in front and end (N+D+N days) to sift through adjustments
def patchData(hist_data:pd.DataFrame(), period:int, debug:bool=0) -> pd.DataFrame():
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

#system constant
eps = np.finfo(float).eps
def chol_corr(Z: np.ndarray) -> np.ndarray:
# compute cholesky decomp of correlation matrix of columns of Z
# attempts to repair non-positive-definite matrices
# http://www.mathworks.com/matlabcentral/answers/6057-repair-non-positive-definite-correlation-matrix
# rank-1 update followed by rescaling to get unit diagonal entries
    R = np.corrcoef(Z)
    posDef = True
    i = 0
    while posDef:
        i += 1
        if i > 1000:
            raise Exception('Reached 1000 iterations, will not continue. Please check time series data')
        
        #will raise LinAlgError if the matrix is not positive definite.
        #https://stackoverflow.com/a/16266736
        try :
            U = np.linalg.cholesky(R)
            posDef = False
        except np.linalg.LinAlgError:
            w,v = np.linalg.eig(R)
            temp = np.amin(v.real)
            k = min(temp,-1*eps)
            R = R - k*np.eye(R.shape[0]) #eye needs 1 dimension only
            R = R/R[0,0] #adjusted for indices
            
            #these doesn'
            #old = R
            #R = get_near_psd(R)
            #R = fix_nonpositive_semidefinite(R)
            #if R.all() == old.all():
            #    print('same')
            posDef = True
        except Exception as err:
            print(err)
            raise Exception('Unable to handle this error')
    return U

def monthly_main( hist_data:pd.DataFrame(), numRealisations:int, numYears:int ) -> {'realisation':pd.DataFrame()} :
    print('Monthly stream flow generation in progress ...')
    # from daily to monthly
    Qh_mon = convert_data_to_monthly(hist_data) 
    Qgen = {}
    for r in range(0,numRealisations):
        Qmon_gen = monthly_gen(Qh_mon, numYears) #formerly Qs
        Qgen[r] = Qmon_gen #raveled in original version
    
    print('done')
    return Qgen #{realisation}{site}{year,month} - years in rows


##-------------------------------------------------------------------------------
##Unverified functions
def monthly_gen(q_historical:pd.DataFrame(), num_years:int, droughtProbab:float=None, nDroughtYears:int=None)-> {'station name':pd.DataFrame()}:
    """
    q_historical is monthly
    """
    qH_stns = list(q_historical.keys()) #keys are station names
    #nPoints = len(qH_stns) #no of donor stations
    nQ_historical = len(q_historical[qH_stns[0]]) #no of years(with months) of data
    #future - check all dfs are of same size in qhistorical
    
    num_years = num_years+1 
    # this adjusts for the new corr technique
    
    if droughtProbab == None and nDroughtYears == None:
        nQ = nQ_historical
    else:
        nDroughtYears = nDroughtYears-1 # (input nDroughtYears=2 to double the frequency, i.e. repmat 1 additional time)
        nQ = nQ_historical + math.floor(droughtProbab*nQ_historical+1)*nDroughtYears
    
    Random_Matrix = np.random.randint(nQ, size=(num_years,12))
    #print(Random_Matrix.shape,'rmsize')
    
    Qs = {}
    for k in q_historical.keys():
        Q_matrix = q_historical[k]
        if  droughtProbab != None and nDroughtYears != None:
            tempDf = Q_matrix.copy()
            tempDf = temp.apply(lambda x: x.sort_values().values)
            # find lowest droughtProbab# of values for each month for drought scenario
            appndDf = tempDf.iloc[0:math.ceil(droughtProbab*nQ_historical),:]
            for i in range(0,nDroughtYears):
                Q_matrix = pd.concat([Q_matrix,appndDf])
                
        logQ = Q_matrix.apply(lambda x: np.log(x))
        #make sure the columns aka months are in sequential order
        logQ = logQ.reindex(sorted(logQ.columns), axis=1)
        
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
        Z = np.array(Z).T #list of lists to array format
        
        Z_vector = np.ravel(Z) #dont transpose here please
        #for across year correlation
        Z_JJ = Z_vector[6:(nQ*12-6)] #:-6 is ideal but drought years are appended, so
        #july of start year to june of end year
        Z_shifted = Z_JJ.reshape(-1,12) #no transpose and check the order too
        #https://bic-berkeley.github.io/psych-214-fall-2016/index_reshape.html
        
        # The correlation matrices should use the historical Z's
        # (the "appendedd years" do not preserve correlation)
        U = chol_corr(Z[:nQ_historical-1,:].T)
        U_shifted = chol_corr(Z_shifted[:nQ_historical-2,:].T)
        
        Qs_uncorr = []
        for i in range(0,12): #python index related
            #print('z RM',Z.shape,Random_Matrix.shape)
            #Qs_uncorr.append(Z[i,Random_Matrix[:,i]])#Z shape changed
            Qs_uncorr.append(Z[Random_Matrix[:,i], i])
        Qs_uncorr = np.array(Qs_uncorr).T #append works in a different way, so .T
        
        Qs_uncorr_vector = np.ravel(Qs_uncorr.T)
        Qs_uncorr_vector_JJ = Qs_uncorr_vector[6:(num_years*12-6)] #index should be correct
        Qs_uncorr_shifted = Qs_uncorr_vector_JJ.reshape(-1,12) #no transpose here too
        
        Qs_corr = np.dot(Qs_uncorr,U)
        Qs_corr_shifted = np.dot(Qs_uncorr_shifted,U_shifted)
        
        Qs_log = np.empty((len( Qs_corr_shifted ),len( Qs_corr_shifted [0])))
        #shifted starts from july, so, pick from 6 month offset to get jan
        Qs_log[:,0:5] = Qs_corr_shifted[:,6:11] #indices adjusted
        #no need of adjustment here
        Qs_log[:,6:11] = Qs_corr[1:num_years,6:11] 
        
        Qsk = np.empty((len( Qs_log ),len( Qs_log[0])))
        for i in range(0,12):
            Qsk[:,i] = np.exp(Qs_log[:,i]*monthly_stdev[i] + monthly_mean[i])
        QskDf = pd.DataFrame(Qsk, columns = Q_matrix.columns)
        Qs[k] = QskDf
        ##print(QskDf.info())

    return Qs

def KNN_identification( Z:{'station name':pd.DataFrame()}, Qtotals:{'offsets':pd.DataFrame()}, year:int, month:int, K:int=None ):
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

    ##Ntotals is calc in two steps here
    wShifts = list(Qtotals.keys()) #offsets
    nSites = list(Qtotals[wShifts[0]].keys()) #stations
    Ntotals = len(wShifts)*len(Qtotals[wShifts[0]][nSites[0]]) #offsets X years
    #this is a short circuit and assumes all are same: valid since internally generated
    
    # Ntotals is the number of historical monthly patterns used for disaggregation.
    # A pattern is a sequence of ndays of daily flows, where ndays is the
    # number of days in the month being disaggregated. Patterns are all
    # historical sequences of length ndays beginning within 7 days before or
    # after the 1st day of the month being disaggregated.

    # nearest neighbors identification
    # only look at neighbors from the same month +/- 7 days
    
    """
    delta = []
    #the sequence doesn't matter, in comparison to original, as they are reordered
    for i in Qtotals.keys(): #shift windows
        for thisStn in Qtotals[i].keys(): #stations
            for thisYear in list((Qtotals[i][thisStn]).index):
                redDf = Qtotals[i][thisStn]
                temp = redDf[redDf.index == thisYear][month].values[0]
                #delta.append(math.pow((temp-Z[thisStn].iloc[year-1,month-1]),2)) #-1 is index context here, month and year loc are assumed
                delta.append(math.pow((temp-Z[thisStn][month][year]),2))
    #print(delta)
    """
    delta = np.zeros((len(wShifts),len(Qtotals[wShifts[0]][nSites[0]]))) #offsets X years
    for thisOffset in Qtotals.keys(): #shift windows
        for thisStn in Qtotals[thisOffset].keys(): #stations
            yrIndex = 0
            for thisYear in list((Qtotals[thisOffset][thisStn]).index):
                redDf = Qtotals[thisOffset][thisStn]
                temp = redDf[redDf.index == thisYear][month].values[0]
                delta[thisOffset,yrIndex] += math.pow((temp-Z[thisStn][month][year]),2)
                yrIndex += 1 #array needs index not a year, not an issue with offset as it is 0~15
    
    Y = pd.Series(np.ravel(np.array(delta)))
    Y_ord = Y.sort_values()
    #K defines the subset representing the sample
    if K == None:
        K = math.floor(math.sqrt(Ntotals))
    KNN_id = Y_ord[0:K] #
    #print(KNN_id)

    # computation of the weights
    f = np.array(range(1,K+1)) #Watch this, dont adjust index to zero, it will yeild infinity
    f1 = 1/f
    #print(f1)
    W = f1 / sum(f1) 
    
    return KNN_id, W, delta
    
def KNN_sampling( KNN_id:pd.Series(), Wcum:np.ndarray, delta:np.ndarray, Qtotals:{'offsets':pd.DataFrame()}, Qdaily:pd.DataFrame(), offsetDaily:pd.DataFrame(), month:int):
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
    r = np.random.rand(1,1)[0,0]
    Wcum =  pd.concat([pd.Series([0]), Wcum], ignore_index=True)
    #print(Wcum)
    
    KNNs = []
    for thisWt,nexWt in zip(Wcum[:-1],Wcum[1:]): 
        #print(r,thisWt,nexWt)
        if (r > thisWt) and (r <= nexWt):
            #one at a time from octave run
            KNNs = Wcum[Wcum == thisWt].index[0] ##something is wrong here
    #print(KNNs)
    thisKNN_id = (KNN_id.copy()).reset_index(drop=True)
    yearID = thisKNN_id[KNNs]
    #print('yid',yearID)
    
    #extract info from Qtotals
    wShifts = list(Qtotals.keys()) #offsets
    nSites = list(Qtotals[wShifts[0]].keys()) #stations
    yList = Qtotals[wShifts[0]][nSites[0]].index.values
    #magic to find the year and offset
    result = np.where(delta == yearID)
    listOfCoordinates = (list(zip(result[0], result[1])))[0]
    thisOffset = wShifts[listOfCoordinates[0]] #shift key to access
    thisYear = yList[listOfCoordinates[1]] #result is 2d tuple now
    #print(thisOffset,thisYear)
        
    # concatenate last 7 days of last year before first 7 days of first year
    # and first 7 days of first year after last 7 days of last year
    nrows = len(Qdaily)
    shifted_Qdaily = (offsetDaily.iloc[thisOffset:thisOffset+nrows,:]).copy() 
    #this offset starts with neg patching as index 0, till pos patching (len + 2X patch)
    shifted_Qdaily.index = Qdaily.index
    #print(shifted_Qdaily.head())
    dailyFlows = shifted_Qdaily[(shifted_Qdaily.index.year == thisYear) & (shifted_Qdaily.index.month == month)].copy()
    
    for thisColumn in dailyFlows.columns:
        dailyFlows[thisColumn] = dailyFlows[thisColumn]/dailyFlows[thisColumn].sum()
    
    #print(dailyFlows.head()) #dataframe of daily flow for corresponding year and month for all stations
    return dailyFlows #df for all stations

def KNNdailyGen(realisation:int, z:{'station':pd.DataFrame()}, Qtotals:{'offsets':pd.DataFrame()}, yrRange:list, monRange:list, hist_data:pd.DataFrame()):
    d = pd.DataFrame()
    startYear = hist_data.index.year.min()
    
    #dailyData offset now, than doing it again and again in loop
    wShifts = list(Qtotals.keys()) #offsets
    patchNDays = int((wShifts[-1])/2) #int casted to remove any oddities
    offsetDaily = patchData(hist_data,patchNDays)
    
    for year in yrRange:
        print('sampling realisation',realisation+1,' for year',year+1, flush = True)
        #sys.stdout.flush()
        for month in monRange:
            [KNN_id, W, delta] = KNN_identification(z, Qtotals, year, month)
            Wcum = pd.Series(np.array(W)).cumsum()
            py = KNN_sampling(KNN_id, Wcum, delta, Qtotals, hist_data, offsetDaily, month)
            #replace the original years in index with predicted year
            py.reset_index(inplace=True)
            #print(py.head())
            py['Date'] = py['Date'].apply(lambda x: x.replace(year = startYear+year))
            py.set_index('Date', inplace=True)
            #temp = py.apply(lambda x: x*(z[x.name].iloc[year,month-1])/(3600*24))
            temp = py.apply(lambda x: x*(z[x.name][month][year])/(3600*24))
            #aka 86400 s/day
            #print('temp',temp.head())
            d = pd.concat([d, temp])
    return [d, realisation]

##---------------------------------------------------------------------------------
def combined_generator(hist_data:pd.DataFrame(), nR:int, nY:int, selectOffsetYears:list = [], selectMonths:list = [], name:'__main__' = None) -> {'realisation':{'station':pd.DataFrame()}} :
    
    #custom date generation or full schedule
    if len(selectOffsetYears) == 0:
        yrRange = range(0,nY)
    else :
        yrRange = selectOffsetYears
    
    if len(selectMonths) == 0:
        monRange = range(1,13)
    else :
        if all(i > 12 for i in selectMonths):
            raise Exception('Month list should not have values between 0,12')
        else: #for the sake of reader
            monRange = selectMonths
        
    hist_data.sort_values(by='Date', inplace=True)
    
    # generation of monthly data via Kirsch et al. (2013):
    # Kirsch, B. R., G. W. Characklis, and H. B. Zeff (2013), 
    # Evaluating the impact of alternative hydro-climate scenarios on transfer 
    # agreements: Practical improvement for generating synthetic streamflows, 
    # Journal of Water Resources Planning and Management, 139(4), 396–406.
    Q_Qg = monthly_main(hist_data, nR, nY ) #dict(realisation) of dicts(sites) of year-mon matrices
    Qh_mon = convert_data_to_monthly(hist_data)
    Qh_stns = list(Qh_mon.keys())
    Nyears = len(Qh_mon[Qh_stns[0]])

    # disaggregation from monthly to daily time step as in Nowak et al. (2010):
    # Nowak, K., Prairie, J., Rajagopalan, B., & Lall, U. (2010). 
    # A nonparametric stochastic approach for multisite disaggregation of 
    # annual to daily streamflow. Water Resources Research, 46(8).

    # Find K-nearest neighbors (KNN) in terms of total monthly flow and 
    # randomly select one for disaggregation. Proportionally scale the flows in
    # the selected neighbor to match the synthetic monthly total. To
    # disaggregate Jan Flows, consider all historical January totals +/- 7
    # days, etc.

    # concatenate last 7 days of last year before first 7 days of first year
    # and first 7 days of first year after last 7 days of last year
    patchNDays = 7
    extra_hist_data = patchData(hist_data,patchNDays)
    
    nrows = len(hist_data)
    Qtotals = {}
    totalShiftDays = patchNDays - (-1*patchNDays) +1 #last value is ignored, so
    for k in range(0,totalShiftDays):
        Qmonthly_shifted = {}
        shifted_hist_data = (extra_hist_data.iloc[k:k+nrows,:]).copy()
        shifted_hist_data.reset_index(inplace=True,drop=True)
        shifted_hist_data.index = hist_data.index
        Qmonthly_shifted = convert_data_to_monthly(shifted_hist_data)
        Qtotals[k] = Qmonthly_shifted #{window}{stations}{year-month}
    
    shelved = False
    try:
        D = shelve.open('predCache', flag='n') #,'c' #alt flag
        shelved = True
    except Exception as er:
        print(er)
        print('using RAM to store data')
        D = {}

    parallel = True
    if name == None:
        parallel = False
        print('serial simulation')
        
    #parallel option
    if parallel:
        import multiprocessing as mp
        pool = mp.Pool(mp.cpu_count()-1) 
        #leave 1 core for the computer to respond
        
        #helper function to put result in correct place for parallel
        def collect_result(results):
            nonlocal D
            [d, r] = results
            d.sort_values(by=['Date'], inplace=True)
            #print(d.head())
            D[str(r)] = d
        if name == "__main__": #this is a hack
            #https://stackoverflow.com/questions/37737085/is-it-possible-to-use-multiprocessing-in-a-module-with-windows
            #https://stackoverflow.com/a/63632161
            for r in range(0,nR):
                z = Q_Qg[r]
                pool.apply_async(KNNdailyGen, args=(r, z, Qtotals, yrRange, monRange, hist_data), callback=collect_result)
    else:
        for r in range(0,nR):
            z = Q_Qg[r]
            [d, _] = KNNdailyGen(r, z, Qtotals, yrRange, monRange, hist_data)
            d.sort_values(by=['Date'], inplace=True)
            #print(d.head())
            D[str(r)] = d
            #print(D[r].head()) #.info())
    
    if parallel:
        pool.close()
        pool.join()  # postpones the execution of next line of code until all processes in the queue are done.
        ##please don't close shelf here, data is required to be returned
        #if shelved:
            #D.close()
    return D, shelved
    


##---------------------------------------------------------------------------------
#lifted functions
##sourced from https://pyportfolioopt.readthedocs.io/en/latest/_modules/pypfopt/risk_models.html
##Martin, R. A., (2021). PyPortfolioOpt: portfolio optimization in Python. Journal of Open Source Software, 6(61), 3066, https://doi.org/10.21105/joss.03066

def fix_nonpositive_semidefinite(matrix, fix_method="spectral"):
    """
    Check if a covariance matrix is positive semidefinite, and if not, fix it
    with the chosen method.

    The ``spectral`` method sets negative eigenvalues to zero then rebuilds the matrix,
    while the ``diag`` method adds a small positive value to the diagonal.

    :param matrix: raw covariance matrix (may not be PSD)
    :type matrix: pd.DataFrame
    :param fix_method: {"spectral", "diag"}, defaults to "spectral"
    :type fix_method: str, optional
    :raises NotImplementedError: if a method is passed that isn't implemented
    :return: positive semidefinite covariance matrix
    :rtype: pd.DataFrame
    """
    
    """
    #know this already
    if _is_positive_semidefinite(matrix):
        return matrix

    warnings.warn(
        "The covariance matrix is non positive semidefinite. Amending eigenvalues."
    )
    """

    # Eigendecomposition
    q, V = np.linalg.eigh(matrix)

    if fix_method == "spectral":
        # Remove negative eigenvalues
        q = np.where(q > 0, q, 0)
        # Reconstruct matrix
        fixed_matrix = V @ np.diag(q) @ V.T
    elif fix_method == "diag":
        min_eig = np.min(q)
        fixed_matrix = matrix - 1.1 * min_eig * np.eye(len(matrix))
    else:
        raise NotImplementedError("Method {} not implemented".format(fix_method))

    if not _is_positive_semidefinite(fixed_matrix):  # pragma: no cover
        warnings.warn(
            "Could not fix matrix. Please try a different risk model.", UserWarning
        )

    # Rebuild labels if provided
    if isinstance(matrix, pd.DataFrame):
        tickers = matrix.index
        return pd.DataFrame(fixed_matrix, index=tickers, columns=tickers)
    else:
        return fixed_matrix
        
#adjust matrix to be positive definite
##https://stackoverflow.com/a/63131309
def get_near_psd(A):
    C = (A + A.T)/2
    eigval, eigvec = np.linalg.eig(C)
    eigval[eigval < 0] = 0

    return eigvec.dot(np.diag(eigval)).dot(eigvec.T)