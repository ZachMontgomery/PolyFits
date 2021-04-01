"""Multivariable polynomial fit using Least Squares Regression.

This module contains functions to calculate and use the polynomial
coefficients for an arbitrary order polynomial curve fit to a dataset with
an arbitrary number of independent variables. Curve fits can be performed
with full control over the polynomial terms and custom weighting of
datapoints.

Methods are taken from:
    Ullah, A. H., Fabijanic, C., Estevadeordal, J., Montgomery, Z. S.,
    Hunsaker, D. F., Staiger, J. M., and Joo, J. J., "Experimental and
    Numerical Evaluation of the Performance of Parabolic Flaps," AIAA
    Aviation 2019 Forum, June 2019,
    https://arc.aiaa.org/doi/abs/10.2514/6.2019-2916
    
    Morelli, E. A., "Global Nonlinear Aerodynamic Modeling using
    Multivariate Orthogonal Functions," Journal of Aircraft, Vol. 32, Issue
    2, 1995, pp. 270-277, https://arc.aiaa.org/doi/abs/10.2514/3.4

Routine Listings
-----------------

multivariablePolynomialFit
    function for calculating a curve fit to data with an arbitrary number of
    independent variables

multivariablePolynomialFunction
    function for calculating a polynomial with an arbitrary number of
    independent variables

multivariableR2
    function calculating the coefficient of determination value, or R^2
    value, of a polynomial fit of an arbitrary number of independent
    variables

multivariableRMS
    function calculating an RMS (root, mean, squared) error and a custom
    RMSN (root, mean, squared, normalized) error where normalized means the
    error is divided by the mean of the absolute value of the dependent
    variables, for a multidimensional polynomial function.

compose_j
    function used by the multivariable series of functions that composes the
    n values of the independent variables into the counter j, this function
    can also be used for the nhat values and i

decompose_j
    function used by the multivariable seris of functions that decomposes
    the counter j into the different n values for the independent variables,
    this function can also be used for i and the nhat values.

calcJ
    Computes the total number of terms in a given multivariable polynomial
    function

kDecompose
    Computes the exponents on the independent variables of the kth term of a
    multivariable polynomial funciton

kCompose
    Computes the k index of the multivariable polynomial term corresponding
    to the exponents of the independent variables.

oneLineProgress
    Class that defines a compact progress bar.

multiSort
    Function that sorts multiple arrays based on a master array, keeping the
    order of data between arrays consistent.

isClose
    Function that determines if two values are 'sufficiently' close together

autoPolyFit
    Performs a multivariate polynomial curve fit to a dataset and
    automatically determines which polynomial terms to use based on a
    balance between the goodness of the fit and a predictve capabilities
    measure that attempts to make the model compact.
    
    Based on the method given by: Morelli, E. A., "Global Nonlinear
    Aerodynamic Modeling using Multivariate Orthogonal Functions," Journal
    of Aircraft, Vol. 32, Issue 2, 1995, pp. 270-277,
    https://arc.aiaa.org/doi/abs/10.2514/3.4

zachsAutoPolyFit
    Under Construction

list2dict
    Function converting the multivariable polynomial function information to
    a dictionary.

dict2list
    Function converting the multivariable polynomial dictionary to the
    polynomial coefficient list, order list, and R&2 value.

polyFit2json
    Function that writes the multivariable polynomial data to a JSON file.
"""
import numpy as np
from datetime import datetime as dt
from datetime import timedelta as td

def multivariablePolynomialFit(Nvec, x, y ,interaction=True, sym=[], sym_same=[], sym_diff=[], zeroConstraints=[], constraints=[], percent=False, weighting=None, calcR2=True, verbose=True):
    """Performs the multivariable polynomial fit using Least Squares Regression.
    
    Parameters
    ----------
    Nvec : list
        List with a length V equal to the number of independent variables.
        The ith element values are integers of polynomial order of the ith
        independent variable.
    
    x : numpy array or nested list
        Numpy array or nested list of size k by V, where k is the total
        number of points in the dataset. The ith column represents the ith
        independent variable values with the rows representing the different
        data points.
    
    y : numpy array or list
        List with a length of k with the dependent variable values
    
    interaction : boolean, optional
        Boolean value with default set to False. This variable determines
        whether or not interaction terms are included in the fit function.
        If set to True, interaction terms up the max order for each
        independent variable are included, i.e. if Nvec = [3,2], then the
        highest interaction term included is x_1^3*x_2^2. Specific
        interaction terms can be omitted using the constraints input.
    
    sym : list, optional
        Optional list that defaults as an empty list. If used, the length
        should be V and each element should contain a boolean, True/False.
        The ith element determines if the ith independent variable is
        symmetric either even or odd, which is determined by the order given
        in Nvec. This will also remove the cooresponding interaction terms
        if they are enabled.
    
    sym_same : list, optional
        Optional list that defaults as an empty list. If used, the entries
        in the list should be tuples with two integers. The integers
        represent the independent variables that the "same" symmetry
        condition will be applied. The "same" symmetry forces all
        polynomial terms between the two independent variables that have the
        same parity on the exponent, i.e. x1^3 * x2^5 both exponents are odd
        so the "same" symmetry would force the poly, coefficient related to
        this term to 0.
    
    sym_diff : list, optional
        Optional list that defaults as an empty list. Same as sym_same,
        except that it implements the "diff" symmetry constraint of opposite
        parity exponent terms being forced to 0.
    
    zeroConstraints : list, optional
        An optional list that defaults as an empty list. Entries in the list
        contain integer tuples of length V. The integer values represent the
        powers of the independent variables whose coefficient will be forced
        to 0 before the best fit calculations are performed, allowing the
        user to omit specific polynomial terms
    
    constraints : list, optional
        An optional list that defaults to an empty list. Entries in the list
        contain tuples of length 2. The first entry is a list of integers
        that represent the powers of the independent variables whose
        coefficient will then be forced to be equal to the second entry in
        the tuple, which should be a float.
    
    percent : boolean, optional
        Boolean value with default set to False. When set to True the least
        squares is performed on the percent error squared. This option
        should not be used if y contains any zero or near zero values, as
        this might cause a divide by zero error.
    
    weighting : callable function, optional
        If given, weighting should be a function that takes as arguments:
        x, y, and p where x and y are the independent and dependent
        variables defined above and p is the index representing a certain
        data point. weighting should return a 'weighting factor' that
        determines how important that datapoint is. Returning a '1' weights
        the datapoint normally, while a '0' will cause the datapoint to be
        ignored.
    
    calcR2 : boolean, optional
        Boolean value with default set to True. Determines if the funciton
        should compute the R^2 value and return it.
    
    verbose : boolean, optional
        Boolean value with default set to True. Determines if the function
        will print relavent data on screen during computation time.
    
    Returns
    -------
    
    list
        List of the polynomial coefficients with a length equal to the
        products of the Nvec elements plus one.
        i.e.: (n_vec0+1)*(n_vec1+1)*...
    
    float, optional
        The coefficient of determination, also referred to as the R^2 value.
        A value of 1 means the polynomial fits the data perfectly.
    """
    ## copy raw data into new arrays
    x = np.copy(x)
    y = np.copy(y)
    
    ## determine number of independent variables
    if type(Nvec) != list:
        Nvec = [Nvec]
    V = len(Nvec)
    
    ## calc number of data points
    k = len(y)
    
    ## check for inconsistencies in dimensions from input variables
    if len(x.shape) == 1:
        x = np.transpose([x])
    if x.shape[1] != V: raise ValueError('Dimensions for V don\'t match between n_vec and x. Length of n_vec and number of columns of x should equal the number of independent variables used, V.')
    if x.shape[0] != k: raise ValueError('Number of rows between x and y don\'t agree! The number of rows should be the total number of points, k, of the dataset.')
    
    ## determine number of polynomial coefficients
    J = calcJ(Nvec)
    
    ## if sym wasn't given, initialize it to False values
    if type(sym) != list: sym = [sym]
    if sym == []:
        sym = [False] * V
    elif len(sym) != V:
        raise ValueError('Length of sym doesn\'t match the number of dimensions, V.')
    
    ## create active list
    ########################################################################
    ## set active to empty list
    active = []
    ## loop through j values
    for j in range(J):
        ## calculate the n values
        n = decompose_j(j, Nvec)
        ## check if n is a constraint then continue on to the next j
        if tuple(n) in zeroConstraints: continue
        ## check if j is an interaction term and if interactions aren't allowed then continue on to the next j
        if sum(n) != max(n) and not interaction: continue
        ## initialize flag variable to false
        flag = False
        ## loop through the sym list to find the symmetry constraints
        for count,symm in enumerate(sym):
            ## check if flag has been tripped, then continue to the next j if it has
            if flag: break
            ## check for a symmetry constraint
            if symm:
                ## check if the order of the count-th independent variable is even
                if Nvec[count]%2 == 0:
                    ## check if the n value of the count-th independent variable is odd
                    if n[count]%2 == 1:
                        flag = True
                ## this else block means the order of the count-th independent variable is odd
                else:
                    ## check if the n value of the count-th independent variable is even
                    if n[count]%2 == 0:
                        flag = True
        ## if the flag has been tripped, skip to the next j value
        if flag: continue
        ## loop through sym_same constraints
        for val in sym_same:
            ## check if the n values from both variables given in val are even, then trip flag
            if n[val[0]]%2 == 0 and n[val[1]]%2 == 0:
                flag = True
            ## check if the n values from both variables given in val are odd, then trip flap
            if n[val[0]]%2 == 1 and n[val[1]]%2 == 1:
                flag = True
        ## loop through sym_diff constraints
        for val in sym_diff:
            ## check if the n values from both variables given in val are even and odd, then trip flag
            if n[val[0]]%2 == 0 and n[val[1]]%2 == 1:
                flag = True
            ## check if the n values from both variables given in val are odd and even, then trip flap
            if n[val[0]]%2 == 1 and n[val[1]]%2 == 0:
                flag = True
        ## if flag hasn't been tripped, append j value onto the active list
        if not flag: active.append(j)
    lenActive = len(active)
    
    ## add weighting factor
    W = np.eye(k)
    if callable(weighting):
        if verbose: prog = oneLineProgress(k, msg='Setting up weighting factors')
        for kk in range(k):
            W[kk,kk] = weighting(x, y, kk) ** 2.
            if verbose: prog.display()
    
    ## set lmda vector
    if percent:
        for kk in range(k):
            W[kk,kk] /= y[kk] ** 2.
    
    ## compute X matrix
    ########################################################################
    ## initialize X matrix to nan's
    X = np.ones((k,lenActive)) * np.nan
    ## set progress bar
    if verbose: prog = oneLineProgress(k*lenActive, msg='PolyFitSetup')
    ## loop thru data points
    for kk in range(k):
        ## loop thru used polynomial coefficients
        for jj,j in enumerate(active):
            ## determine exponents for the polynomial coefficient
            n = decompose_j(j, Nvec)
            ## initialize temp variable to 1
            temp = 1.
            ## loop thru independent variables
            for v in range(V):
                ## multiply on independent variable with the corresponding exponent to temp
                temp *= x[kk,v] ** n[v]
            ## set the temp value to the X matrix
            X[kk,jj] = temp
            if verbose: prog.display()
    
    ## compute A and b
    if verbose: print('Computing the A matrix and b vector')
    A = X.T.dot(W).dot(X)
    b = X.T.dot(W).dot(y)
    
    ## update A and b with the nonzero constraints
    ########################################################################
    for n,val in constraints:
        j = compose_j(n, Nvec)
        # con['{}'.format(j)] = val
        jj = active.index(j)
        for i in range(lenActive):
            if i != jj:
                A[i,jj] = 0.
            else:
                A[i,jj] = 1.
        b[jj] = val
    
    ## solve for the polynomial coefficients
    ########################################################################
    if verbose: print('solving the Aa=b equation')
    a = np.linalg.solve(A,b)
    
    #input the missing 0 coefficients into 'a' so that it can be used with the multidimensional_poly_func
    ########################################################################
    for i in range(J):
        if not i in active:
            a = np.insert(a,i,0.)
            active = np.insert(active,i,0)
    
    a = list(a)
    #calculate R^2 value
    ########################################################################
    if calcR2:
        r = multivariableR2(a, Nvec, x, y, verbose=verbose)
        return a, r
    return a

def multivariablePolynomialFunction(a, Nvec, x):
    """Computes the multivariable polynomial function
    
    Parameters
    ----------
    
    a : list
        List of the polynomial coefficients that has a length equal to
        the products of the Nvec elements plus one.
        i.e.: (Nvec0+1)*(Nvec1+1)*...
    
    Nvec : list
        List with a length V equal to the number of independent
        variables. The ith element values are integers of polynomial
        order of the ith independent variable.
    
    x : list
        List of length V equal to the number of independent variables.
        The ith element represents the ith independent variable values.
    
    Returns
    -------
    
    float
        Value of the multivariable polynomial function for the given
        independent variables
    """
    if type(Nvec) not in (tuple, list, np.ndarray): Nvec = [Nvec]
    if type(x) not in (tuple, list, np.ndarray): x = [x]
    # initialize summation to 0
    f = 0.
    # calculate total number of datapoints
    k = len(a)
    # calculate total number of dimensions
    V = len(x)
    # loop through the datapoints
    for j in range(k):
        # calculate the n values
        n = decompose_j(j, Nvec)
        # initialize the product series variable to 1
        prod = 1.
        # loop through the dimensions
        for v in range(V):
            # multiply onto the product series the term
            prod *= x[v] ** n[v]
        # add onto the summation the proper term
        f += a[j] * prod
    # return the finalized summation value
    return f

def multivariableR2(a, Nvec, x, y, verbose=True):
    """Calculates the R^2 value of a multivariable polynomial fit to a dataset
    
    Parameters
    ----------
    
    a : list
        List of polynomial coefficients
    
    Nvec : list
        List of integers representing the polynomial order of the
        independent variables
    
    x : numpy array or nested list
        Numpy array or nested list of size k by V, where k is the total
        number of points in the dataset. The ith column represents the ith
        independent variable values with the rows representing the different
        data points.
    
    y : numpy array or list
        List with a length of k with the dependent variable values
    
    verbose : boolean, optional
        Determines verbosity of the function. Defaults to True.
    
    Returns
    -------
    
    float
        R^2 value, the coefficient of determination
    """
    # ensure x and y are in the proper format
    x = np.copy(x)
    y = np.copy(y)
    ynew = np.copy([temp for temp in y if temp != None])
    # calculate k
    k = len(ynew)
    # calculate mean y value
    y_ = sum(ynew) / float(k)
    # calculate the SSt value
    SSt = sum( (ynew - y_) ** 2. )
    # initialize the f array
    f = []
    # loop through the datapoints
    if verbose: prog = oneLineProgress(len(y), msg='Determining R^2 of the fit')
    for i in range(len(y)):
        if y[i] == None:
            if verbose: prog.display()
            continue
        # calculate the f value from the polynomial function
        f.append( multivariablePolynomialFunction(a,Nvec,x[i,:]) )
        if verbose: prog.display()
    f = np.copy(f)
    # calculate the SSr term
    SSr = sum( (ynew - f) ** 2. )
    # calculate and return the R^2 value
    return 1. - SSr / SSt

def multivariableRMS(raw_x, raw_y, a, Nvec, verbose=True):
    """Calculates the RMS and RMSN errors of a multivariable polynomial fit to a dataset.
    
    Parameters
    ----------
    
    x : numpy array or nested list
        Numpy array or nested list of size k by V, where k is the total
        number of points in the dataset. The ith column represents the ith
        independent variable values with the rows representing the different
        data points.
    
    y : numpy array or list
        List with a length of k with the dependent variable values
    
    a : list
        List of polynomial coefficients
    
    Nvec : list
        List of integers representing the polynomial order of the
        independent variables
    
    verbose : boolean, optional
        Determines verbosity of the function. Defaults to True.
    
    Returns
    -------
    
    float
        RMS, the root mean squared error
    
    float
        RMSN, the root mean squared of the normalized error
    """
    x = np.copy( raw_x )
    y = np.copy( raw_y )
    avg = np.mean(abs(y))
    k = len(x[:,0])
    func = np.zeros(k)
    e = np.zeros(k)
    e_per = np.zeros(k)
    if verbose: prog = oneLineProgress(k, msg='Determining RMS of the fit')
    for i in range(k):
        func[i] = multivariablePolynomialFunction(a, Nvec, x[i])
        e[i] = (y[i] - func[i]) ** 2.
        e_per[i] = ((y[i] - func[i])/avg) ** 2.
        if verbose: prog.display()
    return np.sqrt(np.mean(e)), np.sqrt(np.mean(e_per))

def compose_j(n, Nvec):
    """Calculates the j index of the multivariable polynomial function.
    
    Parameters
    ----------
    
    n : list
        List of integer values representing the independent variables'
        exponents for the jth term in the multivariable polynomial function.
    
    Nvec : list
        List of integers representing the polynomial order of the
        independent variables.
    
    Returns
    -------
    
    integer
        j index of the multivariable polynomial function
    """
    # calculate V
    V = len(Nvec)
    # initialize j to 0
    j = 0
    # loop through independent variables
    for v in range(1,V+1):
        # initialize product series to 1
        prod = 1
        # loop through w values for product series
        for w in range(v+1,V+1):
            # multiply on the term to the product series
            prod *= Nvec[w-1] + 1
        # add on term onto j
        j += n[v-1] * prod
    return j

def decompose_j(j, Nvec):
    """Calculates the independent variables' exponents corresponding to the jth polynomial coefficient.
    
    Parameters
    ----------
    
    j : integer
        j index representing the jth polynomial term
    
    Nvec : list
        List of integers representing the polynomial order of the
        independent variables.
    
    Returns
    -------
    
    list
        List of integer values representing the independent variables'
        exponents for the jth term in the multidimensional polynomial
        function.
    """
    # calculate V
    V = len(Nvec)
    # initialize n values to nothing
    n = [None]*V
    # loop through the n values that need to be solved, starting at the highest and working down
    for v in range(V,0,-1):
        # initialize the denomenator product series to 1
        denom = 1
        # loop through the w values needed for the product series
        for w in range(v+1,V+1):
            # multiply on the terms for the denomenator product series
            denom *= Nvec[w-1] + 1
        # initialize the summation variable to 0
        summ = 0
        # loop through the u values necessary for the summation
        for u in range(v+1,V+1):
            # initialize the product series variable inside the summation to 1
            prod = 1
            # loop through the s values needed for the product series that is inside the summation
            for s in range(u+1,V+1):
                # multiply on the term for the product series that is inside of the summation
                prod *= Nvec[s-1] + 1
            # add on the needed term to the summation series
            summ += n[u-1] * prod
        # finally calculate the n value cooresponding to this v
        n[v-1] = int(round( ((j-summ)/denom)%(Nvec[v-1]+1) ))
    return n

############################################################################
############################################################################
############################################################################

def calcJ(Nvec):
    '''Computes the total number of terms in a given multivariable polynomial function
    
    Parameters
    ----------
    
    Nvec : list
        List of integers representing the polynomial order of the
        independent variables.
    
    Returns
    -------
    
    integer
        number of terms in the multivariable polynomial function
    '''
    J = 1
    for n in Nvec:
        J *= n + 1
    return J

def kDecompose(k, V):
    '''Computes the exponents on the independent variables of the kth term of a multivariable polynomial funciton
    
    Parameters
    ----------
    
    k : integer
        kth term of a multivariable polynomial function
    
    V : integer
        number of independent variables in the polynomial funciton
    
    Returns
    -------
    
    list
        exponents on the independent variables 
    '''
    t = 1
    ###################################################
    ## find the category
    c = 0
    vals = [0] * V
    while k > t:
        c += 1
        m = [c] * (V-1)
        vals[0] = calcJ(m)
        for j in range(V-1):
            m[j] -= 1
            vals[j+1] = calcJ(m)
        t += sum(vals)
    if c == 0:
        return [0]*V
    ####################################################
    ## find the subcategory
    for sc in range(V-1,-1,-1):
        t -= vals[sc]
        if k > t:
            break
    ####################################################
    ## determine n
    # initialize n
    n = [None]*V
    n[sc] = c
    # create mx based on the sc, then decompose to get m
    mx = [c]*(V-1)
    for i in range(sc):
        mx[i] -= 1
    m = decompose_j(k-t-1, mx)
    # set m values into n and return
    j = -1
    for i in range(V):
        if i != sc:
            j += 1
            n[i] = m[j]
    return n

def kCompose(n):
    '''Computes the k index of the multivariable polynomial term corresponding to the exponents of the independent variables.
    
    Parameters
    ----------
    
    n : list
        exponents of the indepedent variables
    
    Returns
    -------
    
    integer
        k index of the corresponding polynomial term
    '''
    ########################################
    V = len(n)
    if V == 1: return calcJ(n)
    mx = max(n)
    if mx == 0: return 1
    k = 1
    ## calculate lower number sets
    for i in range(1,mx):
        m = [i] * (V-1)
        k += calcJ(m)
        for j in range(V-1):
            m[j] -= 1
            k += calcJ(m)
    ## calculate location in current number set
    for i in range(V):
        M = [mx]*(V-1)
        for j in range(i):
            M[j] -= 1
        if n[i] != mx:
            k += calcJ(M)
        else:
            m = [n[j] for j in range(V) if j != i]
            k += compose_j(m, M) + 1
            return k
    raise ValueError('Unable to compose n into k: current k value {}'.format(k))

class oneLineProgress():
    '''A class defining a compact progress bar
    
    Parameters
    ----------
    
    Total : integer
        number of iterations for the progress event
    
    msg : string, optional
        Message displayed for the progress bar. Defaults to no message.
    
    showETR : boolean, optional
        Determines whether to display the estimated time remaining. Defaults to True.
    
    Returns
    -------
    
    oneLineProgress
        a newly initialized compact progress bar
    '''
    
    def __init__(self, Total, msg='', showETR=True):
        self.total = Total
        self.msg = msg
        self.count = 0
        self.showETR = showETR
        self.start = dt.now()
        self.rollTimer = dt.now()
        self.rollCount = -1
        self.rollDelta = 0.2
        self.display()
    
    def increment(self):
        '''Increments the counter variable for the progress bar.'''
        self.count += 1
    
    def decrement(self):
        '''Decrements the counter variable for the progress bar.'''
        self.count -= 1
    
    def __str__(self):
        '''Under Construction'''
        pass
    
    def __len__(self):
        '''Under Construction'''
        l = len(str(self))
        self.decrement()
        return l
    
    def Set(self, count):
        '''Set the value of the couter variable for the progress bar.'''
        self.count = count
    
    def display(self):
        '''Increments the counter and updates the progress bar if needed.'''
        rolling = '-\\|/'
        rollDelta = (dt.now()-self.rollTimer).total_seconds()
        
        p2s = False
        if rollDelta >= self.rollDelta or self.rollCount == -1:
            p2s = True
            self.rollTimer = dt.now()
            self.rollCount += 1
            if self.rollCount >= len(rolling):
                self.rollCount = 0
        
        perc = self.count / self.total * 100.
        self.increment()
        
        if not p2s and perc < 100.: return
        
        s = '\r' + ' '*(len(self.msg)+50) + '\r'
        s += self.msg + ' '*4
        
        # j = 0
        for i in range(10):
            if perc >= i*10:
                j = i
        
        if perc < 100.:
            s += u'\u039e'*j + rolling[self.rollCount] + '-'*(9-j)
        else:
            s += u'\u039e'*10
        
        # for i in range(1,11):
            # if i*10 <= perc:
                # s += u'\u039e'
            # else:
                # s += '-'
        s += ' '*4 + '{:7.3f}%'.format(perc)
        if not self.showETR:
            if perc >= 100.: s += '\n'
            print(s, end='')
            return
        
        if perc <= 0:
            etr = '-:--:--.------'
            s += ' '*4 + 'ETR = {}'.format(etr)
        elif perc >= 100.:
            s += ' '*4 + 'Run Time {}'.format(dt.now()-self.start) + '\n'
        else:
            time = (dt.now()-self.start).total_seconds()
            etr = td(seconds=time / perc * 100. - time)
            s += ' '*4 + 'ETR = {}'.format(etr)
        print(s, end='')
        return

def multiSort(v, *W, ascend=True, verbose=True, msg='Sorting the arrays'):
    '''Sorts multiple arrays based on a master array.
    
    Parameters
    ----------
    
    v : list or numpy array
        Master array to be sorted
    
    W : list(s) or numpy array(s) or combinations, optional
        Additional lists of the same length as the master that will be
        sorted to match the order of the master
    
    ascend : boolean, optional
        Determines if the master array should be sorted in ascending order.
        Defaults to True.
    
    verbose : boolean, optional
        Determines if a progress bar will update the status of the sorting.
        Defaults to True.
    
    msg : string, optional
        Message to be included with the progress bar. Defaults to 'Sorting
        the arrays'.
    '''
    k = len(v)
    for w in W:
        if len(w) != k: raise ValueError('All arrays need to be the same length in multiSort')
    c = []
    if verbose: prog = oneLineProgress(sum([i for i in range(k)])+len(W), msg=msg)
    for m in range(k):
        for j in range(k-1,m,-1):
            i = j-1
            
            if (ascend and v[j] < v[i]) or (not ascend and v[j] > v[i]):
                c.append(j)
                temp = v[j]
                v[j] = v[i]
                v[i] = temp
            if verbose: prog.display()
    
    for w in W:
        for j in c:
            i = j-1
            temp = w[j]
            w[j] = w[i]
            w[i] = temp
        if verbose: prog.display()

def isClose(x, y, tol=1.e-12):
    '''Determines if two values are sufficiently close together.
    
    Parameters
    ----------
    
    x : float, int, or other '<' '>' comparable object
        First value
    
    y : same as x
        Second value. The order of x and y doesn't matter.
    
    tol : same as x, optional
        Tolerance value that determines what is 'sufficiently close'.
        Defaults to 1.e-12.
    
    Returns
    -------
    
    boolean
        True if the two values are close, otherwise False.
    '''
    return y-tol <= x and x <= y+tol


def autoPolyFit(X, y, MaxOrder=12, tol=1.e-12, sigma=None, sigmaMultiplier=1., verbose=True):
    '''Automatic Multivariable Polynomial Curve Fit
    
    Performs a multivariate polynomial curve fit to a dataset and
    automatically determines which polynomial terms to use based on a
    balance between the goodness of the fit and a predictve capabilities
    measure that attempts to make the model compact.
    
    Based on the method given by: Morelli, E. A., "Global Nonlinear
    Aerodynamic Modeling using Multivariate Orthogonal Functions," Journal
    of Aircraft, Vol. 32, Issue 2, 1995, pp. 270-277,
    https://arc.aiaa.org/doi/abs/10.2514/3.4
    
    Parameters
    ----------
    X : numpy array
        Array of shape (N,m). X consists of all the independent variables in
        the dataset. N is the number of data points in the set and m is the
        number of independent variables
    
    y : list or numpy array
        Array with length N. y is the dependent variable values
        cooresponding to the independent variables in X
    
    MaxOrder : integer, optional
        Gives the max order of polynomial for any one of the independent
        varialbes to try. Defaults to 12
    
    tol : float, optional
        Gives the cut-off value for any polynomial coefficient to not be
        included in the final results. If a coefficient has an absolute
        value below tol, it won't be included. Defaults to 1e-12
    
    sigma : float, optional
        Value used to determine the trade off between how good of a fit to
        perform and how many terms to keep. Defaults to None, which causes
        the function to calculate sigma automatically using the mean squared
        of the difference of the independent variable values with respect to
        the mean independent variable value of the dataset
    
    sigmaMultiplier : float, optional
        Term multiplied onto sigma to change it's value. Allows using a
        multiple of the automatically determined sigma value. Defaults to 1.
    
    verbose : boolean, optional
        Determines the verbosity of the function. Defaults to True.
    
    Returns
    -------
    list
        A list of the polynomial coefficients.
    
    list
        A list of the max polynomial orders for each independent variable.
        The length of this list is therefore m. This list is comparable to
        the 'Nvec' object used throughout this module
    
    float
        The coefficient of determination, R^2 value, representing the
        goodness of the fit
    '''
    ## number of independent variables
    m = X.shape[1]
    ## max range of polynomials to try
    Nvec = tuple([MaxOrder]*m)
    ## number of datapoints
    N = len(y)
    ## number of p functions
    K = kCompose(Nvec)
    ###################################################################
    ##           determine the orthogonal p functions
    if verbose: prog = oneLineProgress(K-1, msg='Determining the orthogonal p functions')
    ## initialize the P matrix
    P = np.zeros((N, K))
    P[:,0] = 1.
    ## loop thru k values
    for k in range(2,K+1):
        ## determine the set associated with k
        n = decompose_j(k-1, Nvec)
        ## find pkhat and mu
        mu = None
        for i in range(m):
            nhat = n[:]
            if nhat[i] > 0:
                nhat[i] -= 1
                khat = compose_j(nhat, Nvec) + 1
                if khat < k:
                    mu = i
                    break
        if mu == None: raise ValueError('Unable to find previous khat set')
        pkhat = P[:,khat-1]
        xmu = X[:,mu]
        ## calculate pk first term
        temp = xmu * pkhat
        phik = sum(n)
        pk = temp[:]
        ## loop thru summation in eq 18
        for j in range(1,k):
            ## check if value is needed
            phij = sum(decompose_j(j-1, Nvec))
            if phik - phij <= 2:
                ## calculate gamma
                pj = P[:,j-1]
                zzz = sum([abs(i) for i in pj])
                if zzz == 0.: print(zzz)
                gamma = np.dot(pj, temp) / np.dot(pj, pj)
                pk -= gamma * pj
        ## add pk to P
        P[:,k-1] = pk[:]
        if verbose: prog.display()
    #################################################################
    ##              sort the p functions by effectiveness
    order = [i for i in range(K)]
    ranks = [None] * K
    for i in range(K):
        pj = P[:,i]
        pjdot = np.dot(pj, pj)
        ajhat = np.dot(pj, y) / pjdot
        ranks[i] = ajhat ** 2. * pjdot
    multiSort(ranks, order, ascend=False, msg='Sorting the p functions by effectivenss', verbose=verbose)
    Pordered = np.zeros((N,K))
    for i,o in enumerate(order):
        Pordered[:,i] = P[:,o]
    P = Pordered[:,:]
    ###################################################################
    ##          determine how many of the orthogonal p functions to use
    if verbose: prog = oneLineProgress(K, msg='Determining number of p functions to use')
    PSEold = None
    foundMin = False
    if sigma == None:
        yavg = sum(y) / N
        sigma = sum([(i - yavg)**2. for i in y]) / N
    sigma *= sigmaMultiplier
    
    for n in range(1,K+1):
        Phat = P[:,:n]
        ahat = np.matmul(np.matmul(np.linalg.inv(np.matmul(Phat.transpose(), Phat)), Phat.transpose()), y)
        yhat = np.matmul(Phat, ahat)
        MSE = np.dot(y - yhat, y - yhat) / N
        PSEnew = MSE + sigma * n / N
        if verbose: prog.display()
        if PSEold == None or PSEnew <= PSEold:
            PSEold = PSEnew
        else:
            foundMin = True
            nn = n-1
            P = Phat[:,:nn]
            order = order[:nn]
            if verbose: 
                prog.Set(K)
                prog.display()
            break
    if not foundMin:
        # raise ValueError('Unable to find minimum PSE')
        print('========== Unable to find minimum PSE ==========')
        nn = K
    ###################################################################
    ##              final coefficients and polynomial size
    if verbose: prog = oneLineProgress(4+nn, msg='Determining final coefficients and polynomial size')
    b = np.zeros((nn,nn))
    for k in range(1,nn+1):
        j = k - 1
        pj = P[:,j]
        w = np.ones((N,k))
        for i in range(k):
            n = decompose_j(order[i], Nvec)
            for ii,e in enumerate(n):
                w[:,i] *= X[:,ii] ** e
        vals = np.matmul(np.matmul(np.linalg.inv(np.matmul(w.transpose(), w)), w.transpose()), pj)
        b[j,:k] = vals[:]
        if verbose: prog.display()
    
    A = [np.dot(P[:,i],y)/np.dot(P[:,i],P[:,i]) for i in range(nn)]
    if verbose: prog.display()
    
    c = [np.dot(A,b[:,i]) for i in range(nn)]
    js = [decompose_j(order[i], Nvec) for i in range(nn)]
    if verbose: prog.display()
    
    js = [js[i] for i in range(nn) if not isClose(c[i], 0., tol=tol)]
    c = [i for i in c if not isClose(i, 0., tol=tol)]
    if verbose: prog.display()
    
    nvec = [None]*m
    for i in range(m):
        nvec[i] = max([j[i] for j in js])
    JJ = 1
    for n in nvec:
        JJ *= n+1
    
    a = [0.] * JJ
    for j in range(JJ):
        n = decompose_j(j, nvec)
        for i in range(len(c)):
            if n == js[i]:
                a[j] = c[i]
    if verbose: prog.display()
    
    return a, nvec, multivariableR2(a, nvec, X, y, verbose=verbose)


def zachsAutoPolyFit(X, y, ranges, MaxOrder=12, tol=1.e-2, verbose=True):
    '''Under Construction'''
    V = X.shape[1]
    nvec = [MaxOrder]*V
    Cons, cons = [], []
    c = 1
    while c == 1 or len(cons) > 0:
        c += 1
        
        a = multivariablePolynomialFit(nvec, X, y, interaction=True, zeroConstraints=Cons, calcR2=False, verbose=verbose)
        cons = []
        J = calcJ(nvec)
        for j in range(J):
            n = decompose_j(j,nvec)
            
            temp = 1
            for i in range(V):
                temp *= abs(ranges[i])**n[i]
            temp = tol/temp #* 10.
            
            if isClose(a[j], 0., tol=temp):
                if a[j] != 0:
                    cons.append(tuple(n))
                    print(n)
                # a[j] = 0.
        # print(cons)
        # print('Count = {}'.format(len(cons)))
        Cons += cons
    
    js = [decompose_j(j,nvec) for j in range(J) if a[j] != 0.]
    Nvec = [None]*V
    for i in range(V):
        Nvec[i] = max([j[i] for j in js])
    JJ = calcJ(Nvec)
    A = [0.]*JJ
    for j in range(J):
        if a[j] != 0.:
            n = decompose_j(j,nvec)
            i = compose_j(n,Nvec)
            A[i] = a[j]
    
    # a = A[:]
    # nvec = Nvec[:]
    
    
    R2 = multivariableR2(A, Nvec, X, y, verbose=verbose)
    return A, Nvec, R2

def list2dict(a, Nvec, R2):
    '''Returns a dictionary containing the multivariable polynomial function information.'''
    return {'Nvec':Nvec, 'Coefficients':list(a), 'R^2':R2}

def dict2list(data):
    '''Returns the multivariable polynomial coefficients list, order list, and R^2 value from a multivariable polynomial dictionary.'''
    return data['Coefficients'], data['Nvec'], data['R^2']

def polyFit2json(filename, a, nvec, r2):
    '''Saves the multivariable polynomial funtion information into a JSON file.'''
    import json
    d = list2dict(a, nvec, r2)
    f = open(filename, 'w')
    json.dump(d, f, indent=4)
    f.close()
    return d

