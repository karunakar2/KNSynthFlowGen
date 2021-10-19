#import sys

import math
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
    Nsites = len(Qh_mon)

    Qgen = {}
    for r in range(0,nR):
        Qgen[r] = {}
        Qmon_gen = monthly_gen(Qh_mon, nY) #formerly Qs
        Qgen[r] = Qmon_gen
    
    #print(Qgen,'Qgen')
    ##print(len(Qgen),Qgen[1].shape)
    return Qgen #{realisation}{site}{year,month} - years in rows

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
        
        Z_vector = np.ravel(Z.T)
        Z_JulyToJune = Z_vector[6:-6] #zero index so starts at 6, uptill last 6 months
        #july of start year to june of end year
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
            Qs_uncorr.append(Z[i,Random_Matrix[:,i]])
        Qs_uncorr = np.array(Qs_uncorr).T
        #print('c',Qs_uncorr.shape)
        Qs_corr = np.dot(Qs_uncorr,U)
        
        Qs_uncorr_vector = np.ravel(Qs_uncorr)
        Qs_uncorr_vector_JJ = Qs_uncorr_vector[6:-6]
        Qs_uncorr_shifted = Qs_uncorr_vector_JJ.reshape(-1,12)
        Qs_corr_shifted = np.dot(Qs_uncorr_shifted,U_shifted)

        Qs_log = np.empty((len( Qs_corr_shifted ),len( Qs_corr_shifted [0])))
        #print(Qs_log)
        Qs_log[:,0:5] = Qs_corr_shifted[:,6:11] #indices adjusted
        Qs_log[:,6:11] = Qs_corr[:-1,6:11] #not sure should we exclude last or first?
        #print(Qs_log.shape)

        Qsk = np.empty((len( Qs_log ),len( Qs_log[0])))
        for i in range(0,12):
            Qsk[:,i] = np.exp(Qs_log[:,i]*monthly_stdev[i] + monthly_mean[i])
        QskDf = pd.DataFrame(Qsk, columns = Q_matrix.columns)
        Qs[k] = QskDf
        ##print(QskDf.info())

    return Qs

###---------------------yet to finish this part -----------------------------------

def KNN_identification( Z, Qtotals, year, month, K=None ):
    #year and month are being used in index context here
    #quick shortcut, as they both are numeric
    ##year -= 1
    ##month -= 1
    
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
    wShifts = list(Qtotals.keys())
    nSites = list(Qtotals[wShifts[0]].keys())
    Ntotals = len(wShifts)*len(Qtotals[wShifts[0]][nSites[0]]) #this is a short circuit
    
    # Ntotals is the number of historical monthly patterns used for disaggregation.
    # A pattern is a sequence of ndays of daily flows, where ndays is the
    # number of days in the month being disaggregated. Patterns are all
    # historical sequences of length ndays beginning within 7 days before or
    # after the 1st day of the month being disaggregated.
    if K == None:
        K = math.floor(math.sqrt(Ntotals))

    # nearest neighbors identification
    # only look at neighbors from the same month +/- 7 days
    delta = []
    for i in Qtotals.keys(): #shift windows
        for j in Qtotals[i].keys(): #stations
            for thisYear in list((Qtotals[i][j]).index):
                temp = Qtotals[i][j]
                temp = temp[temp.index == thisYear][month].values[0]
                delta.append(math.pow((temp-Z[j].iloc[year-1,month-1]),2)) #-1 is index context here
    #print(delta)
    Y = pd.Series(np.array(delta))
    Y_ord = Y.sort_values()
    KNN_id = Y_ord[0:K] #
    #print(KNN_id)

    # computation of the weights
    f = np.array(range(1,K+1)) #Watch this, dont adjust index to zero, it will yeild infinity
    f1 = 1/f
    #print(f1)
    W = f1 / sum(f1) 

    return KNN_id, W, Y_ord
    
def KNN_sampling( KNN_id, Wcum, Y_ord, Qdaily, month, windowSize):
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
            #KNNs.append(i) #is this really supposed to be more than 1? or a single value
            #one at a time from octave run
            KNNs = Wcum[Wcum == thisWt].index[0] ##something is wrong here
    #print(KNNs)
    thisKNN_id = (KNN_id.copy()).reset_index(drop=True)
    yearID = thisKNN_id[KNNs]
    #print('yid',yearID)
    
    #magic to find the year and k
    KNNmat = np.array(Y_ord).reshape(windowSize,len(Qdaily.columns),-1)
    result = np.where(KNNmat == yearID)
    listOfCoordinates = list(zip(result[0], result[1], result[2])) #3d array
    k = listOfCoordinates[0][0] #use the first hit
    if k >= windowSize:
        raise Exception('The number of days shifted is more than set window size ', windowSize)
    year = listOfCoordinates[0][2]
    #see KNN_identification() for the format and size

    # concatenate last 7 days of last year before first 7 days of first year
    # and first 7 days of first year after last 7 days of last year
    nrows = len(Qdaily)
    Qdaily = patchData(Qdaily,k)
      
    shifted_Qdaily = (Qdaily.iloc[k:k+nrows,:]).copy()
    #print(2,len(Qdaily),len(shifted_Qdaily))
    shifted_Qdaily.index = Qdaily.index
    #print(shifted_Qdaily.head())
    thisYear = int(shifted_Qdaily.index.year.min()) + year
    #print(thisYear)
    dailyFlows = shifted_Qdaily[(shifted_Qdaily.index.year == thisYear) & (shifted_Qdaily.index.month == month)].copy()
    
    for thisColumn in dailyFlows.columns:
        dailyFlows[thisColumn] = dailyFlows[thisColumn]/dailyFlows[thisColumn].sum()
    
    #print(dailyFlows.head()) #dataframe of daily flow for corresponding year and month for all stations
    return dailyFlows #df for all stations

def KNNdailyGen(realisation, z, Qtotals, yrRange, monRange, hist_data, totalShiftDays):
    d = pd.DataFrame()
    for year in yrRange:
        print('sampling realisation',realisation+1,' for year',year+1, flush = True)
        #sys.stdout.flush()
        for month in monRange:
            [KNN_id, W, Y_ord] = KNN_identification(z, Qtotals, year, month)
            Wcum = pd.Series(np.array(W)).cumsum()
            py = KNN_sampling(KNN_id, Wcum, Y_ord, hist_data, month, totalShiftDays)
            temp = py.apply(lambda x: x*z[x.name].iloc[year,month-1]) #/Dt) 
            # division with Dt is reducing the values by several orders
            d = pd.concat([d, temp])
    return [d, realisation]

##---------------------------------------------------------------------------------
#DaysPerMonth = [31,28,31,30,31,30,31,31,30,31,30,31]
def combined_generator(hist_data, nR, nY, selectYears = [], selectMonths = [], name = None) :
    
    #custom date generation or full schedule
    if len(selectYears) == 0:
        yrRange = range(0,nY)
    else :
        yrRange = selectYears
        
    if len(selectMonths) == 0:
        monRange = range(1,13)
    else :
        monRange = selectMonths
    
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
    totalShiftDays = patchNDays - (-1*patchNDays)
    for k in range(0,totalShiftDays):
        Qmonthly_shifted = {}
        shifted_hist_data = (extra_hist_data.iloc[k:k+nrows,:]).copy()
        shifted_hist_data.reset_index(inplace=True,drop=True)
        shifted_hist_data.index = hist_data.index
        Qmonthly_shifted = convert_data_to_monthly(shifted_hist_data)
        Qtotals[k] = Qmonthly_shifted #{window}{stations}{year-month}
    
    Dt = 3600*24
    D = {}
    """
    parallel = False
    #print(name)
    if name != None :
        if name == "__main__" or name == "__mp_main__": #former is original call, later is mp call
            parallel = True
            print('running parallel simulation')
        else:
            print('parallel execution is attempted but didnot work because the file from where this function is called, should be named "main.py"')
    if parallel == False:
         print('serial simulation')
    """
    parallel = True
    if name == None:
        parallel = False
        print('serial simulatiom')
        
    #parallel option
    if parallel:
        import multiprocessing as mp
        pool = mp.Pool(mp.cpu_count()-1) 
        #leave 1 core for the computer to respond
        
        #helper function to put result in correct place for parallel
        def collect_result(results):
            global D
            [d, r] = results
            D[r] = d
        if name == "__main__": #this is a hack
            #https://stackoverflow.com/questions/37737085/is-it-possible-to-use-multiprocessing-in-a-module-with-windows
            #https://stackoverflow.com/a/63632161
            for r in range(0,nR):
                z = Q_Qg[r]
                pool.apply_async(KNNdailyGen, args=(r, z, Qtotals, yrRange, monRange, hist_data, totalShiftDays), callback=collect_result)
    else:
        for r in range(0,nR):
            z = Q_Qg[r]
            [D[r], _] = KNNdailyGen(r, z, Qtotals, yrRange, monRange, hist_data, totalShiftDays)
            #print(D[r].head()) #.info())
    
    if parallel:
        pool.close()
        pool.join()  # postpones the execution of next line of code until all processes in the queue are done.
    return D
    


##---------------------------------------------------------------------------------
