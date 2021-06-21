import ZachsModules as zm
import numpy as np

class polyFit():
    
    def __init__(self, Nvec, **kw):
        
        if type(Nvec) != list: raise TypeError()
        for N in Nvec: if type(N) != int: raise TypeError()
        self.Nvec               = Nvec[:]
        self.V                  = len(Nvec)
        self.J                  = self.calcJ()
        
        self.coef               = kw.get('coef', [None]*self.J)
        if type(self.coef) != list: raise TypeError()
        for a in self.coef: if type(a) not in (float, None): raise TypeError()
        
        self.interaction        = kw.get('interaction', True)
        self.sym                = kw.get('sym', [False]*self.V)
        if type(self.sym) != list: self.sym = [self.sym]
        if len(self.sym) != self.V: raise ValueError()
        self.symSame            = kw.get('symSame', [])
        self.symDiff            = kw.get('symDiff', [])
        self.zeroConstraints    = kw.get('zeroContstraints', [])
        self.constraints        = kw.get('constraints', [])
        
        self.Jtilde             = kw.get('Jtilde', None)
        
        self.percent            = kw.get('percent', False)
        self.weighting          = kw.get('weighting', None)
        
        self.R2                 = kw.get('R2', None)
        self.RMS                = kw.get('RMS', None)
        self.RMSN               = kw.get('RMSN', None)
        self.Syx                = kw.get('Syx', None)
        
        self.ybar               = kw.get('ybar', None)
        self.St                 = kw.get('St', None)
        self.Sr                 = kw.get('Sr', None)
    
    def fit(self, x, y, verbose=True):
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
        
        ## calc number of data points
        k = len(y)
        
        ## check for inconsistencies in dimensions from input variables
        if len(x.shape) == 1: x = np.transpose([x])
        if x.shape[1] != self.V: raise ValueError('Dimensions for V don\'t match between n_vec and x. Length of n_vec and number of columns of x should equal the number of independent variables used, V.')
        if x.shape[0] != k: raise ValueError('Number of rows between x and y don\'t agree! The number of rows should be the total number of points, k, of the dataset.')
        
        ## create active list
        ########################################################################
        ## set active to empty list
        active = []
        ## loop through j values
        for j in range(self.J):
            ## calculate the n values
            n = self.decompose_j(j)
            ## check if n is a constraint then continue on to the next j
            if tuple(n) in self.zeroConstraints: continue
            ## check if j is an interaction term and if interactions aren't allowed then continue on to the next j
            if not self.interaction and sum(n) != max(n): continue
            ## initialize flag variable to false
            flag = False
            ## loop through the sym list to find the symmetry constraints
            for count,symm in enumerate(self.sym):
                ## check if flag has been tripped, then continue to the next j if it has
                if flag: break
                ## check for a symmetry constraint
                if symm:
                    ## check if the order of the count-th independent variable is even
                    if self.Nvec[count]%2 == 0:
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
            for val in self.symSame:
                ## check if the n values from both variables given in val are even, then trip flag
                if n[val[0]]%2 == 0 and n[val[1]]%2 == 0:
                    flag = True
                ## check if the n values from both variables given in val are odd, then trip flap
                if n[val[0]]%2 == 1 and n[val[1]]%2 == 1:
                    flag = True
            ## loop through sym_diff constraints
            for val in self.symDiff:
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
        if callable(self.weighting):
            if verbose: prog = zm.io.oneLineProgress(k, msg='Setting up weighting factors')
            for kk in range(k):
                W[kk,kk] = weighting(x, y, kk) ** 2.
                if verbose: prog.display()
        
        ## set lmda vector
        if self.percent:
            for kk in range(k):
                W[kk,kk] /= y[kk] ** 2.
        
        ## compute X matrix
        ########################################################################
        ## initialize X matrix to nan's
        X = np.ones((k,lenActive)) * np.nan
        ## set progress bar
        if verbose: prog = zm.io.oneLineProgress(k*lenActive, msg='PolyFitSetup')
        ## loop thru data points
        for kk in range(k):
            ## loop thru used polynomial coefficients
            for jj,j in enumerate(active):
                ## determine exponents for the polynomial coefficient
                n = self.decompose_j(j)
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
        if verbose: zm.io.oneLineText('Computing the A matrix and b vector')
        A = X.T.dot(W).dot(X)
        b = X.T.dot(W).dot(y)
        
        ## update A and b with the nonzero constraints
        ########################################################################
        for n,val in self.constraints:
            j = self.compose_j(n)
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
        if verbose: zm.io.oneLineText('solving the Aa=b equation')
        a = np.linalg.solve(A,b)
        
        #input the missing 0 coefficients into 'a' so that it can be used with the multidimensional_poly_func
        ########################################################################
        for i in range(self.J):
            if not i in active:
                a = np.insert(a,i,0.)
                active = np.insert(active,i,0)
        
        self.coef = list(a)
    
    def goodnessParams(self, x, y, verbose=True):
        # ensure x and y are in the proper format
        x = np.copy(x)
        y = np.copy(y)
        ynew = np.copy([temp for temp in y if temp != None])
        # calculate k
        k = len(ynew)
        # calculate mean y value
        self.ybar = sum(ynew) / float(k)
        # calculate the SSt value
        self.St = sum( (ynew - y_) ** 2. )
        # initialize the f array
        f = []
        # loop through the datapoints
        if verbose: prog = zm.io.oneLineProgress(len(y), msg='Determining R^2 of the fit')
        for i in range(len(y)):
            if y[i] == None:
                if verbose: prog.display()
                continue
            # calculate the f value from the polynomial function
            f.append( self.evaluate(x[i,:]) )
            if verbose: prog.display()
        f = np.copy(f)
        # calculate the SSr term
        self.Sr = sum( (ynew - f) ** 2. )
        # calculate the R^2 value
        self.R2 = 1. - SSr / SSt
        
        self.Syx = np.sqrt(self.Sr / (k - self.Jtilde))
    

class database():
    
    def __init__(self, ):
        
        self.X = 
        self.Y = 
        
        self.namesX = 
        self.namesY = 
        
        self.V = 
        self.k = 
        
        

