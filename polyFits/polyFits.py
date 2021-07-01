import numpy as np
import ZachsModules as zm
from multiprocessing import Pool, cpu_count

class database():
    
    def __init__(self, x, y, namesX=(), namesY=()):
        
        ## copy input data as numpy array
        self.x = np.copy(x)
        self.y = np.copy(y)
        
        ## check for 1D data and adjust as necessary
        if len(self.x.shape) == 1: self.x = np.transpose([self.x])
        if len(self.y.shape) == 1: self.y = np.transpose([self.y])
        
        ## calc number of data points
        self.numPoints = self.x.shape[0]
        if self.y.shape[0] != self.numPoints: raise ValueError('Number of data points do not match between x and y data inputs with {} and {}'.format(self.x.shape[0], self.y.shape[0]))
        
        self.numIndVar = self.x.shape[1]
        self.numDepVar = self.y.shape[1]
        
        if namesX == ():
            self.namesX = tuple(['x{}'.format(i+1) for i in range(self.numIndVar)])
        else:
            self.namesX = namesX
        if namesY == ():
            self.namesY = tuple(['y{}'.format(i+1) for i in range(self.numDepVar)])
        else:
            self.namesY = namesY

class polyFit():
    
    def __init__(self, db, **kw):
        
        self.db = db
        
        self.Nvec               = kw.get('Nvec', [0]*self.db.numIndVar)
        
        if type(self.Nvec) != list: raise TypeError()
        for N in self.Nvec:
            if type(N) != int: raise TypeError()
        
        self.numCoef            = self.calcNumCoef()
        
        self.coef               = kw.get('coef', [None]*self.numCoef)
        
        if type(self.coef) != list: raise TypeError()
        for a in self.coef:
            if type(a) not in (float, type(None)): raise TypeError('a is type {} but should be either float or None'.format(type(a)))
        
        self.interaction        = kw.get('interaction', True)
        self.sym                = kw.get('sym', [False]*self.db.numIndVar)
        if type(self.sym) != list: self.sym = [self.sym]
        if len(self.sym) != self.db.numIndVar: raise ValueError()
        
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
        
        self.mp                 = kw.get('mp', None)
    
    def createX(self, args):
        kk, jj, j = args
        ## determine exponents for the polynomial coefficient
        n = self.decompose_j(j)
        ## initialize temp variable to 1
        temp = 1.
        ## loop thru independent variables
        for v in range(self.db.numIndVar):
            ## multiply on independent variable with the corresponding exponent to temp
            val = self.db.x[kk,v]
            for _ in range(n[v]): temp *= val
            # temp *= x[kk,v] ** n[v]
        return kk, jj, temp
    
    def fit(self, verbose=True, computeGoodnessParams=True):
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
        
        k = self.db.numPoints
        
        ## create active list
        ########################################################################
        ## set active to empty list
        active = []
        ## loop through j values
        for j in range(self.numCoef):
            ## calculate the n values
            n = self.decompose_j(j)
            ## check if n is a zero constraint then continue on to the next j
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
                if flag: break
                ## check if the n values from both variables given in val are even, then trip flag
                if n[val[0]]%2 == 0 and n[val[1]]%2 == 0:
                    flag = True
                ## check if the n values from both variables given in val are odd, then trip flap
                if n[val[0]]%2 == 1 and n[val[1]]%2 == 1:
                    flag = True
            if flag: continue
            ## loop through sym_diff constraints
            for val in self.symDiff:
                if flag: break
                ## check if the n values from both variables given in val are even and odd, then trip flag
                if n[val[0]]%2 == 0 and n[val[1]]%2 == 1:
                    flag = True
                ## check if the n values from both variables given in val are odd and even, then trip flap
                if n[val[0]]%2 == 1 and n[val[1]]%2 == 0:
                    flag = True
            ## if flag hasn't been tripped, append j value onto the active list
            if not flag: active.append(j)
        lenActive = len(active)
        
        ## compute X matrix
        ########################################################################
        ## initialize X matrix to nan's
        X = np.ones((k,lenActive)) * np.nan
        ## set progress bar
        if verbose: prog = zm.io.oneLineProgress(k*lenActive, msg='PolyFit Setup: Computing the Basis Functions')
        
        if self.mp == 1:
            ## loop thru data points
            for kk in range(k):
                ## loop thru used polynomial coefficients
                for jj,j in enumerate(active):
                    X[kk,jj] = self.createX((kk, jj, j))[-1]
                    if verbose: prog.display()
        else:
            
            it = [None]*(k*lenActive)
            cnt = -1
            ## loop thru data points
            for kk in range(k):
                ## loop thru used polynomial coefficients
                for jj,j in enumerate(active):
                    cnt += 1
                    it[cnt] = (kk, jj, j)
            
            if self.mp == None:
                cpus = cpu_count()
            else:
                cpus = self.mp
            
            with Pool(cpus) as pool:
                chu = lenActive*k//cpus//20
                if chu > 8000: chu = 8000
                if chu < 1: chu = 1
                for vals in pool.imap_unordered(self.createX, it, chunksize=chu):
                    kk, jj, val = vals
                    X[kk, jj] = val
                    if verbose: prog.display()
        
        
        ## add weighting factor
        if callable(self.weighting) or self.percent:
            
            try:
                W = np.eye(k)
                
                if verbose:
                    if callable(self.weighting) and self.percent:
                        pLen = k*2
                    else:
                        pLen = k
                    prog = zm.io.oneLineProgress(k, msg='Setting up weighting factors')
                
                if callable(self.weighting):
                    for kk in range(k):
                        W[kk,kk] = weighting(self.db, kk) ** 2.
                        if verbose: prog.display()
                
                ## set lmda vector
                if self.percent:
                    for kk in range(k):
                        W[kk,kk] /= self.db.y[kk] ** 2.
                        if verbose: prog.display()
                
                ## compute A and b
                if verbose: zm.io.oneLineText('Computing the A matrix and b vector')
                
                A = X.T.dot(W).dot(X)
                b = X.T.dot(W).dot(self.db.y)
                
                del W
                
            except MemoryError:
                
                if verbose:
                    zm.io.oneLineText('Insufficient memory available for matrix implementation of weighting.')
                    prog = zm.io.oneLineProgress(k, msg='Setting up weighting factors manually')
                
                A = np.zeros((lenActive,lenActive))
                b = np.zeros(lenActive)
                
                # W = [1.] * k
                # if callable(self.weighting): W = [self.weighting(self.db, kk) for kk in range(k)]
                # W = [
                
                Xt = X.T.copy()
                
                for kk in range(k):
                    w = 1.
                    if callable(self.weighting): w *= weighting(self.db, kk) ** 2.
                    if self.percent: w /= self.db.y[kk] **2.
                    Xt[:,kk] *= w
                    prog.display()
                
                zm.io.oneLineText('Computing the A matrix and b vector')
                
                A = Xt.dot(X)
                b = Xt.dot(self.db.y)
        else:
            if verbose: zm.io.oneLineText('Computing the A matrix and b vector')
            A = X.T.dot(X)
            b = X.T.dot(self.db.y)
        
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
        
        self.Jtilde = lenActive - len(self.constraints)
        
        ## solve for the polynomial coefficients
        ########################################################################
        if verbose: zm.io.oneLineText('solving the Aa=b equation')
        a = np.linalg.solve(A,b)
        
        #input the missing 0 coefficients into 'a' so that it can be used with the multidimensional_poly_func
        ########################################################################
        for i in range(self.numCoef):
            if not i in active:
                a = np.insert(a,i,0.)
                active = np.insert(active,i,0)
        
        self.coef = list(a)
        
        if computeGoodnessParams:
            self.goodnessParams(verbose=verbose)
    
    def goodnessParams(self, verbose=True):
        ## check for holes
        ynew = np.copy([temp for temp in self.db.y[:,0] if temp != None])
        # calculate k
        k = len(ynew)
        # calculate mean y value
        self.ybar = sum(ynew) / float(k)
        # calculate the SSt value
        self.St = float(sum( (ynew - self.ybar) ** 2. ))
        # initialize the f array
        f = np.zeros(k)
        # loop through the datapoints
        if verbose: prog = zm.io.oneLineProgress(self.db.numPoints, msg='Evaluating Fit Parameters')
        
        if self.mp == 1:
            cnt = 0
            for i in range(self.db.numPoints):
                if self.db.y[i] == None:
                    continue
                ## calculate the f value from the polynomial function
                f[cnt] = self.evaluate(self.db.x[i,:])
                cnt += 1
                if verbose: prog.display()
        else:
            it = []
            for i in range(self.db.numPoints):
                if self.db.y[i] == None: continue
                it.append(self.db.x[i,:])
            
            if self.mp == None:
                cpus = cpu_count()
            else:
                cpus = self.mp
            
            with Pool(cpus) as pool:
                chu = k // cpus // 20
                if chu > 8000: chu = 8000
                if chu < 1: chu = 1
                for i, val in enumerate(pool.imap(self.evaluate, it, chunksize=chu)):
                    f[i] = val
                    if verbose: prog.display()
        
        # calculate the SSr term
        self.Sr = float(sum( (ynew - f) ** 2. ))
        # calculate the R^2 value
        self.R2 = 1. - self.Sr / self.St
        
        if k - self.Jtilde != 0:
            self.Syx = np.sqrt(self.Sr / (k - self.Jtilde))
        else:
            self.Syx = float('nan')
        
        avg = np.mean(abs(ynew))
        e = np.zeros(k)
        e_per = np.zeros(k)
        # for i in range(k):
            # e[i] = (y[i] - f[i]) ** 2
        
        # self.RMS  = np.sqrt(np.mean((ynew - f) ** 2.))
        self.RMS = np.sqrt(self.Sr / k)
        self.RMSN = np.sqrt(np.mean(((ynew - f)/avg) ** 2.))
    
    def evaluate(self, x):
        if type(x) not in (tuple, list, np.ndarray): x = [x]
        if len(x) != self.db.numIndVar: raise ValueError()
        # initialize summation to 0
        f = 0.
        # loop through the coefficients
        for j in range(self.numCoef):
            # calculate the n values
            n = self.decompose_j(j)
            # initialize the product series variable to 1
            prod = 1.
            # loop through the dimensions
            for v in range(self.db.numIndVar):
                # multiply onto the product series the term
                prod *= x[v] ** n[v]
            # add onto the summation the proper term
            f += self.coef[j] * prod
        # return the finalized summation value
        return f
    
    def compose_j(self, n):
        # initialize j to 0
        j = 0
        # loop through independent variables
        for v in range(1,self.db.numIndVar+1):
            # initialize product series to 1
            prod = 1
            # loop through w values for product series
            for w in range(v+1,self.db.numIndVar+1):
                # multiply on the term to the product series
                prod *= self.Nvec[w-1] + 1
            # add on term onto j
            j += n[v-1] * prod
        return j
    
    def decompose_j(self, j):
        # initialize n values to nothing
        n = [None]*self.db.numIndVar
        # loop through the n values that need to be solved, starting at the highest and working down
        for v in range(self.db.numIndVar,0,-1):
            # initialize the denomenator product series to 1
            denom = 1
            # loop through the w values needed for the product series
            for w in range(v+1,self.db.numIndVar+1):
                # multiply on the terms for the denomenator product series
                denom *= self.Nvec[w-1] + 1
            # initialize the summation variable to 0
            summ = 0
            # loop through the u values necessary for the summation
            for u in range(v+1,self.db.numIndVar+1):
                # initialize the product series variable inside the summation to 1
                prod = 1
                # loop through the s values needed for the product series that is inside the summation
                for s in range(u+1,self.db.numIndVar+1):
                    # multiply on the term for the product series that is inside of the summation
                    prod *= self.Nvec[s-1] + 1
                # add on the needed term to the summation series
                summ += n[u-1] * prod
            # finally calculate the n value cooresponding to this v
            n[v-1] = int(round( ((j-summ)/denom)%(self.Nvec[v-1]+1) ))
        return n
    
    def calcNumCoef(self):
        J = 1
        for n in self.Nvec:
            J *= n + 1
        return J
    
    def __str__(self):
        return zm.io.text([ 'R^2 = {}'.format(self.R2),
                            'RMS = {}'.format(self.RMS),
                            'RMSN = {}'.format(self.RMSN),
                            'Syx = {}'.format(self.Syx)], p2s=False)

