import numpy as np
import ZachsModules as zm
from multiprocessing import Pool, cpu_count
import os
import shutil
import json

class database():
    
    def __init__(self, x, y, name='database1', namesX=(), namesY=()):
        
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
        
        self.name = name

class polyFit():
    
    def __init__(self, db, kw, mpFits=1, verbose=True):
        '''Performs polynomial fits to a data set.
        
        Parameters
        ----------
        db : instance of the database class
            containes the independent and dependent variable values of the
            data set
        
        kw : list
            list of length db.numDepVar equal to the number of dependent
            variables in the data set. The i-th entry of the list contains a
            dictionary of keyword arguments for the fit of the ith dependent
            variable.
            
            There are two types of fits that can be performed, each of which
            contain their own keyword arguments. The first is an automatic
            polynomial fit and the second is a manual polynomial fit.
            
            The automatic polynomial fit is based on:
            Morelli, E. A., "Global Nonlinear Aerodynamic Modeling using
            Multivariate Orthogonal Functions," Journal of Aircraft, Vol 32,
            Issue 2, 1995, pp. 270-277,
            https://arc.aiaa.org/doi/abs/10.2514/3.4.
            This option performs a multivariable polynomial fit to the data
            set and automatically determines which polynomial terms to use
            based on a balance between the goodness of the fit and a
            predictive capabilities measure that attempts to make the model
            compact.
            
            The kw arguments used for this type of fit are:
                maxOrder : int, optional
                    gives the max order of polynomial for any one of the
                    independent variables to try, defaults to 6
                
                tol : float, optional
                    gives the cut-off value for any polynomial coefficient
                    to not be included in the final results. If a coef. has
                    an absolute value below tol, it will not be included in
                    the final polynomial coefficients, defaults to 1e-12
                
                sigma : float, optional
                    value used to determine the trade off between how good
                    of a fit to perform and how many terms to keep. Defaults
                    to None, which causes the function to calculate sigma
                    automatically using the mean squared of the difference
                    of the independent variable values with respect to the
                    mean independent variable value of the data set
                
                sigmaMultiplier : float, optional
                    term multiplied onto sigma to change it's value. Allows
                    using a multiple of the automatically determined sigma
                    value, defaults to 1.
                
                verbose : boolean, optional
                    determines the verbosity of the function, defaults to 
                    True.
            
            The manual polynomial fit is based on:
            Ullah, A. H., Fabijanic, C., Estevadeordal, J., Montgomery, Z.
            S., Hunsaker, D. F., Staiger, J. M., and Joo, J. J.,
            "Experimental and Numerical Evaluation of the Performance of
            Parabolic Flaps," AIAA Aviation 2019 Forum, June 2019,
            https://arc.aiaa.org/doi/abs/10.2514/6.2019-2916.
            This option performs a multivariable polynomial fit to the data
            set for a given set of polynomial orders of each of the
            independent variables and cen be performed with full control of
            the polynomial fit by allowing: various symmetry constraints,
            forcing given coefficients to zero or a set value, and custom
            weighting of the data points.
            
            The kw arguments for this type of fit are:
                Nvec : list
                    list of ints that is length db.numIndVar equal to the
                    number of independent variables, the ith element gives
                    the polynomial order of the ith independent variable
                
                interaction : boolean, optional
                    determines whether or not interaction terms are included
                    in the fit function., if set to True, interation terms
                    up to the max order for each independent variable are
                    included, i.e. if Nvec = [3,2], then the highest
                    interaction term included is x_1^3*x_2^2. Specific
                    interaction terms (and other terms) can be omitted using
                    other options. Defaults to True
                
                sym : list, optional
                    list of length db.numIndVar, the number of independent
                    variables where each element is a boolean. The ith entry
                    determines whether the ith independent variable has a
                    symmetry constraint applied. The type of symmetry, even
                    or odd, is determined by the max polynomial order for
                    that variable given in Nvec. Defaults to a list of False
                
                crossSymEven : list, optional
                    list of tuples containing two integers. The two integers            double check this
                    represent the independent variables that the even cross
                    symmetry is applied, which is an even symmetry that acts
                    along a diagonal among the two independent variables 
                    listed. The even cross symmetry forces all polynomial
                    terms between the two independent variables that have
                    the same parity on the exponents, i.e. x1^3 * x2^5 both
                    exponents are odd so the even cross symmetry would force
                    the poly coefficient related to this term to 0. Defaults
                    to an empty list
                
                crossSymOdd : list, optional
                    list of tuples containing two integers. The two integers            double check this
                    represent the independent variables that the odd cross
                    symmetry is applied, which is an odd symmetry that acts
                    along a diagonal among the two independent variables 
                    listed. The odd cross symmetry forces all polynomial
                    terms between the two independent variables that have
                    different parity on the exponents, i.e. x1^4 * x2^5 have
                    exponents are odd and even, so the odd cross symmetry
                    would force the poly coefficient related to this term to
                    0. Defaults to an empty list
                
                zeroConstraints : list, optional
                    entries in the list are tuples of ints of length equal
                    to the number of independent variables. The ints
                    represent the exponents on the independent variables of
                    the polynomial term whose coefficient will be forced to
                    0 before the best fit calculations are performed,
                    allowing the user to omit specific polynomial terms.
                    Defaults to an empty list
                
                constraints : list, optional
                    list of tuples each of length 2. The first entry in the
                    tuple is a list of integers that represent the exponents
                    on the independent variables of the polynomial term,
                    whose coefficient will be forced to the second entry in
                    the tuple, which should be float. This constraint is
                    applied before the best fit is calculated, ensuring that
                    the fit equations take into account the constraint.
                    Defaults to an empty list
                
                percent : boolean, optional
                    when set to True, it performs the least squares fit on
                    the percent error instead of just the error. This option
                    should not be used if any of the dependent varialbe
                    values are 0 or near zero as this might cause division
                    by zero errors. This effectively causes the method to
                    weight data points of lower magnitude as more important
                    than data ponts of higher magnitude. Defaults to False
                
                weighting : callable function, optional
                    If given, weighting should be a function that takes as
                    arguments db, iy, and ip where db is the database object
                    iy is the index representing the independent variable, 
                    and ip is the index representing a certain data point.
                    weighting should return a 'weighting factor' that
                    determines how important that datapoint is. Returning a
                    '1' weights the datapoint normally (assuming the average
                    of all the weighting factors returned is 1), while a '0'
                    will cause the datapoint to be ignored.
                
                mp : int, optional
                    determines the number of cpus to use to perform the poly
                    fit. Setting to 1 causes the function to not use any 
                    multiprocessing. Defaults to the number of cpus on your
                    machine
                
                verbose : boolean, optional
                    determines the verbosity of the function, defaults to 
                    True.
            
            mpFits : int, optional
                determines the number cpus to allow the fits to be computed
                simultaneously using multiprocessing. If set to 1,
                multiprocessing will not be used and fits of the dependent
                variables will be performed one at a time. If set to a value
                other than 1, multiprocessing will be used and 'mp' from the
                manual fit methods will be overridden to be 1 so as to not
                conflict with this multiprocessing.
                
                Note, than any manual fits that have all the exact same 'kw'
                arguments can and will be run simultaneously without the
                need of multiprocessing.
                
                defaults to 1
                
                if set to 0, the number of cpus on your machine will be used
                
                disables verbosity for each individual fit
            
            verbose : boolean, optional
                determines the verbosity of the major steps of the data base
                curve fitting process.
        '''
        
        ## copy in database
        self.db = db
        
        ## create array for polynomial values corresponding to the db values
        self.f = np.zeros((self.db.numPoints, self.db.numDepVar))
        
        ## copy in fit arguments
        self.kw = kw[:]
        
        ## determine which fits will use auto
        self.auto = [False] * self.db.numDepVar
        for i in range(self.db.numDepVar):
            if self.kw[i].get('Nvec', None) == None: self.auto[i] = True
        
        ## initialize the global Nvec
        self.Nvec    = [None] * self.db.numDepVar
        
        ## initialize the global number of coefficients
        self.numCoef = np.array( [0]*self.db.numDepVar )
        
        ## initialize the global coefficient list
        self.coef    = [None] * self.db.numDepVar
        
        ## initialize other global goodness measurement variables
        self.Jtilde  = np.array( [0]*self.db.numDepVar )
        self.R2      = np.zeros( self.db.numDepVar )
        self.RMS     = np.zeros( self.db.numDepVar )
        self.RMSN    = np.zeros( self.db.numDepVar )
        self.Syx     = np.zeros( self.db.numDepVar )
        self.ybar    = np.zeros( self.db.numDepVar )
        self.St      = np.zeros( self.db.numDepVar )
        self.Sr      = np.zeros( self.db.numDepVar )
        
        ###################################
        ## check for duplicate manual fits
        ###################################
        ## initialize an empty list to track duplicates
        self.duplicateManFits = []
        ## initialize a list of the unused dep var indices for the manual fits
        unusedManFits = [i for i in range(self.db.numDepVar) if not self.auto[i]]
        ## loop while there are still untracked manual fits
        while len(unusedManFits) > 0:
            ## initialize a duplicates list with the next untracked man fit
            duplicates = [unusedManFits.pop(0)]
            ## loop thru the remaining untracked man fits
            for i in range(len(unusedManFits)-1,-1,-1):
                ## check if the untracked man fit is the same as man fit of focus
                if self.kw[unusedManFits[i]] == self.kw[duplicates[0]] and not self.kw[duplicates[0]].get('percent', False):
                    ## track the match
                    duplicates.append( unusedManFits.pop(i) )
            ## sort the duplicates
            zm.nm.zSort(duplicates, verbose=False)
            ## add the duplicates to the global variable
            self.duplicateManFits.append( tuple(duplicates) )
        # ## loop thru the fits and print the duplicates
        # if verbose:
            # for dupFit in self.duplicateManFits:
                # k = len(dupFit)
                # if k > 1:
                    # zm.io.text('Found duplicate fit(s) for variables:', *[' {}'.format(self.db.namesY[i]) for i in dupFit])
        
        ## initialize a list of the auto fit dep var indices
        autoFits = [i for i in range(self.db.numDepVar) if self.auto[i]]
        
        ###################
        ## Perform the fits
        ###################
        ## check if perfoming fits individually
        if mpFits == 1:
            ## perform the manual fits
            for dupFit in self.duplicateManFits:
                if verbose: zm.io.text('Performing manual fit(s) for:', *['{}'.format(self.db.namesY[i]) for i in dupFit])
                self.manFit(dupFit)
            ## peform the auto fits
            for i in autoFits:
                if verbose: zm.io.text('Performing auto fit for {}'.format(self.db.namesY[i]))
                self.autoFit(i)
        else: ## perfoming the fits simultanuously with multiprocessing
            ## disable the multiprocessing option for all manual fits
            for i in range(self.db.numDepVar):
                if not self.auto[i]: self.kw[i]['mp'] = 1
            ## initialize iterable for the multiprocessing
            it = []
            ## add on the args for the manual fits
            for dupFit in self.duplicateManFits:
                it.append( dupFit )
            ## add on the args for the auto fits
            for i in autoFits:
                it.append( i )
            ## check if using the same numer of cpus as on the current machine
            if mpFits == 0: mpFits = cpu_count()
            ## initialize progress bar
            if verbose: prog = zm.io.Progress(len(it), title='Performing fits for {}'.format(self.db.name))
            ## perform the fits
            with Pool(mpFits) as pool:
                for _ in pool.imap_unordered(self.whichFit, it):
                    if verbose: prog.display()
        
        ## display the results
        if verbose: print(self)
    
    ####################################################################
    ####################################################################
    ####################################################################
    
    def whichFit(self, iy):
        if type(iy) != int:
            self.manFit(iy)
        else:
            self.autoFit(iy)
        return
    
    def createX(self, args):
        kk, jj, iy, j = args
        ## determine exponents for the polynomial coefficient
        n = self.decompose_j(j, self.Nvec[iy])
        ## initialize temp variable to 1
        temp = 1.
        ## loop thru independent variables
        for v in range(self.db.numIndVar):
            ## multiply on independent variable with the corresponding exponent to temp
            val = self.db.x[kk,v]
            for _ in range(n[v]): temp *= val
            # temp *= x[kk,v] ** n[v]
        return kk, jj, temp
    
    def manFit(self, iy):
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
        
        ## initialize useful variables
        k = self.db.numPoints
        z = iy[0]
        
        ## unpack kw
        Nvec            = self.kw[z]['Nvec']
        interaction     = self.kw[z].get('interaction', True)
        sym             = self.kw[z].get('sym', [False]*self.db.numIndVar)
        crossSymEven    = self.kw[z].get('crossSymEven', [])
        crossSymOdd     = self.kw[z].get('crossSymOdd', [])
        zeroConstraints = self.kw[z].get('zeroConstraints', [])
        constraints     = self.kw[z].get('constraints', [])
        percent         = self.kw[z].get('percent', False)
        weighting       = self.kw[z].get('weighting', None)
        mp              = self.kw[z].get('mp', cpu_count())
        verbose         = self.kw[z].get('verbose', True)
        
        ## set global variables numCoef and Nvec
        self.numCoef[iy,] = J = self.calcNumCoef(Nvec)
        for i in iy: self.Nvec[i] = Nvec[:]
        
        ## create active list
        ########################################################################
        ## set active to empty list
        active = []
        ## loop through j values
        for j in range(J):
            ## calculate the n values
            n = self.decompose_j(j, Nvec)
            ## check if n is a zero constraint then continue on to the next j
            if tuple(n) in zeroConstraints: continue
            ## check if j is an interaction term and if interactions aren't allowed then continue on to the next j
            if not interaction and sum(n) != max(n): continue
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
            for val in crossSymEven:
                if flag: break
                ## check if the n values from both variables given in val are even, then trip flag
                if n[val[0]]%2 == 0 and n[val[1]]%2 == 0:
                    flag = True
                ## check if the n values from both variables given in val are odd, then trip flap
                if n[val[0]]%2 == 1 and n[val[1]]%2 == 1:
                    flag = True
            if flag: continue
            ## loop through sym_diff constraints
            for val in crossSymOdd:
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
        
        if mp == 1:
            ## loop thru data points
            for kk in range(k):
                ## loop thru used polynomial coefficients
                for jj,j in enumerate(active):
                    X[kk,jj] = self.createX((kk, jj, z, j))[-1]
                    if verbose: prog.display()
        else:
            
            it = [None]*(k*lenActive)
            cnt = -1
            ## loop thru data points
            for kk in range(k):
                ## loop thru used polynomial coefficients
                for jj,j in enumerate(active):
                    cnt += 1
                    it[cnt] = (kk, jj, z, j)
            
            if mp == 0:
                cpus = cpu_count()
            else:
                cpus = mp
            
            with Pool(cpus) as pool:
                chu = lenActive*k//cpus//20
                if chu > 8000: chu = 8000
                if chu < 1: chu = 1
                for vals in pool.imap_unordered(self.createX, it, chunksize=chu):
                    kk, jj, val = vals
                    X[kk, jj] = val
                    if verbose: prog.display()
        
        
        ## add weighting factor
        if callable(weighting) or percent:
            
            try:
                W = np.eye(k)
                
                if verbose:
                    if callable(weighting) and percent:
                        pLen = k*2
                    else:
                        pLen = k
                    prog = zm.io.oneLineProgress(pLen, msg='Setting up weighting factors')
                
                if callable(weighting):
                    for kk in range(k):
                        W[kk,kk] = weighting(self.db, z, kk) ** 2.
                        if verbose: prog.display()
                
                ## set lmda vector
                if percent:
                    for kk in range(k):
                        W[kk,kk] /= self.db.y[kk,z] ** 2.
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
                
                Xt = X.T.copy()
                
                for kk in range(k):
                    w = 1.
                    if callable(weighting): w *= weighting(self.db, z, kk) ** 2.
                    if self.percent: w /= self.db.y[kk, z] **2.
                    Xt[:,kk] *= w
                    prog.display()
                
                zm.io.oneLineText('Computing the A matrix and b vector')
                
                A = Xt.dot(X)
                b = Xt.dot(self.db.y)
                
                del Xt
        else:
            if verbose: zm.io.oneLineText('Computing the A matrix and b vector')
            A = X.T.dot(X)
            b = X.T.dot(self.db.y)
        
        ## update A and b with the nonzero constraints
        ########################################################################
        for n,val in constraints:
            j = self.compose_j(n, Nvec)
            jj = active.index(j)
            for i in range(lenActive):
                if i != jj:
                    A[i,jj] = 0.
                else:
                    A[i,jj] = 1.
            b[jj] = val
        
        self.Jtilde[iy,] = lenActive - len(constraints)
        
        ## solve for the polynomial coefficients
        ########################################################################
        if verbose: zm.io.oneLineText('solving the Aa=b equation')
        a = np.linalg.solve(A,b)
        
        ## extract coefficinets
        for i in iy:
            self.coef[i] = a[:,i]
        
        #input the missing 0 coefficients into 'a' so that it can be used with the multidimensional_poly_func
        ########################################################################
        for j in range(J):
            if not j in active:
                for i in iy:
                    self.coef[i] = np.insert(self.coef[i],j,0.)
                active = np.insert(active,j,0)
        
        for i in iy:
            self.coef[i] = list(self.coef[i])
        
        self.goodnessParams(mp, iy, verbose=verbose)
    
    def goodnessParams(self, mp, iy, verbose=True):
        for z in iy:
            ## check for holes
            ynew = np.copy([temp for temp in self.db.y[:,z] if temp != None])
            ## calculate k
            k = len(ynew)
            ## calculate mean y value
            self.ybar[z] = sum(ynew) / float(k)
            ## calculate the SSt value
            self.St[z] = float(sum( (ynew - self.ybar[z]) ** 2. ))
            
            ## loop through the datapoints
            if verbose: prog = zm.io.oneLineProgress(k, msg='Evaluating Fit Parameters for {}'.format(self.db.namesY[z]))
            
            if mp == 1:
                cnt = 0
                for i in range(self.db.numPoints):
                    if self.db.y[i,z] == None:
                        continue
                    ## calculate the f value from the polynomial function
                    self.f[i,z] = self.evaluate(z, self.db.x[i,:])
                    cnt += 1
                    if verbose: prog.display()
            else:
                it = []
                for i in range(self.db.numPoints):
                    if self.db.y[i,z] == None: continue
                    it.append((z,self.db.x[i,:]))
                
                if mp == 0:
                    cpus = cpu_count()
                else:
                    cpus = mp
                
                with Pool(cpus) as pool:
                    chu = k // cpus // 20
                    if chu > 8000: chu = 8000
                    if chu < 1: chu = 1
                    for i, val in enumerate(pool.starmap(self.evaluate, it, chunksize=chu)):
                        self.f[i,z] = val
                        if verbose: prog.display()
            
            # calculate the SSr term
            self.Sr[z] = float(sum( (ynew - self.f[:,z]) ** 2. ))
            # calculate the R^2 value
            self.R2[z] = 1. - self.Sr[z] / self.St[z]
            
            if k - self.Jtilde[z] != 0:
                self.Syx[z] = np.sqrt(self.Sr[z] / (k - self.Jtilde[z]))
            else:
                self.Syx[z] = float('nan')
            
            avg = np.mean(abs(ynew))
            e = np.zeros(k)
            e_per = np.zeros(k)
            # for i in range(k):
                # e[i] = (y[i] - f[i]) ** 2
            
            # self.RMS  = np.sqrt(np.mean((ynew - f) ** 2.))
            self.RMS[z]  = np.sqrt(self.Sr[z] / k)
            self.RMSN[z] = np.sqrt(np.mean(((ynew - self.f[:,z])/avg) ** 2.))
    
    def evaluate(self, z, x):
        if type(x) not in (tuple, list, np.ndarray): x = [x]
        if len(x) != self.db.numIndVar: raise ValueError()
        # initialize summation to 0
        f = 0.
        # loop through the coefficients
        for j in range(self.numCoef[z]):
            # calculate the n values
            n = self.decompose_j(j, self.Nvec[z])
            # initialize the product series variable to 1
            prod = 1.
            # loop through the dimensions
            for v in range(self.db.numIndVar):
                # multiply onto the product series the term
                prod *= x[v] ** n[v]
            # add onto the summation the proper term
            f += self.coef[z][j] * prod
        # return the finalized summation value
        return f
    
    '''
    @staticmethod
    def compose_j(iy, n):
        # initialize j to 0
        j = 0
        # loop through independent variables
        for v in range(1,self.db.numIndVar+1):
            # initialize product series to 1
            prod = 1
            # loop through w values for product series
            for w in range(v+1,self.db.numIndVar+1):
                # multiply on the term to the product series
                prod *= self.Nvec[iy][w-1] + 1
            # add on term onto j
            j += n[v-1] * prod
        return j
    
    @staticmethod
    def decompose_j(self, iy, j):
        # initialize n values to nothing
        n = [None]*self.db.numIndVar
        # loop through the n values that need to be solved, starting at the highest and working down
        for v in range(self.db.numIndVar,0,-1):
            # initialize the denomenator product series to 1
            denom = 1
            # loop through the w values needed for the product series
            for w in range(v+1,self.db.numIndVar+1):
                # multiply on the terms for the denomenator product series
                denom *= self.Nvec[iy][w-1] + 1
            # initialize the summation variable to 0
            summ = 0
            # loop through the u values necessary for the summation
            for u in range(v+1,self.db.numIndVar+1):
                # initialize the product series variable inside the summation to 1
                prod = 1
                # loop through the s values needed for the product series that is inside the summation
                for s in range(u+1,self.db.numIndVar+1):
                    # multiply on the term for the product series that is inside of the summation
                    prod *= self.Nvec[iy][s-1] + 1
                # add on the needed term to the summation series
                summ += n[u-1] * prod
            # finally calculate the n value cooresponding to this v
            n[v-1] = int(round( ((j-summ)/denom)%(self.Nvec[iy][v-1]+1) ))
        return n
    '''
    
    @staticmethod
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
    
    @staticmethod
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
    
    @staticmethod
    def calcNumCoef(Nvec):
        J = 1
        for n in Nvec:
            J *= n + 1
        return J
    
    def __printFit__(self, z):
        msg =           [   'R^2 = {}'.format(self.R2[z]),
                            'RMS = {}'.format(self.RMS[z]),
                            'RMSN = {}'.format(self.RMSN[z]),
                            'Syx = {}'.format(self.Syx[z])]
        if self.auto[z]: msg.insert(0, 'Nvec = {}'.format(self.Nvec[z]))
        return zm.io.text(msg, p2s=False, title=self.db.namesY[z])
    
    def __str__(self):
        s = ''
        for z in range(self.db.numDepVar): s += self.__printFit__(z)
        return s
    
    ####################################################################
    ####################################################################
    ####################################################################
    
    def kDecompose(self, k):
        V = self.db.numIndVar
        t = 1
        ###################################################
        ## find the category
        c = 0
        vals = [0] * V
        while k > t:
            c += 1
            m = [c] * (V-1)
            vals[0] = self.calcNumCoef(m)
            for j in range(V-1):
                m[j] -= 1
                vals[j+1] = self.calcNumCoef(m)
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
        m = self.decompose_j(k-t-1, mx)
        # set m values into n and return
        j = -1
        for i in range(V):
            if i != sc:
                j += 1
                n[i] = m[j]
        return n
    
    @staticmethod
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
        if V == 1: return polyFit.calcNumCoef(n)
        mx = max(n)
        if mx == 0: return 1
        k = 1
        ## calculate lower number sets
        for i in range(1,mx):
            m = [i] * (V-1)
            k += polyFit.calcNumCoef(m)
            for j in range(V-1):
                m[j] -= 1
                k += polyFit.calcNumCoef(m)
        ## calculate location in current number set
        for i in range(V):
            M = [mx]*(V-1)
            for j in range(i):
                M[j] -= 1
            if n[i] != mx:
                k += polyFit.calcNumCoef(M)
            else:
                m = [n[j] for j in range(V) if j != i]
                k += polyFit.compose_j(m, M) + 1
                return k
        raise ValueError('Unable to compose n into k: current k value {}'.format(k))
    
    def autoFit(self, z): #X, y, MaxOrder=12, tol=1.e-12, sigma=None, sigmaMultiplier=1., verbose=True):
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
        
        maxOrder        = self.kw[z].get('maxOrder', 6)
        tol             = self.kw[z].get('tol', 1e-12)
        sigma           = self.kw[z].get('sigma', None)
        sigmaMultiplier = self.kw[z].get('sigmaMultiplier', 1.)
        verbose         = self.kw[z].get('verbose', True)
        mp              = self.kw[z].get('mp', 1)
        
        m = self.db.numIndVar
        ## max range of polynomials to try
        Nvec = tuple([maxOrder]*m)
        ## number of datapoints
        N = self.db.numPoints
        ## number of p functions
        K = self.kCompose(Nvec)
        ###################################################################
        ##           determine the orthogonal p functions
        if verbose: prog = zm.io.oneLineProgress(K-1, msg='Determining the orthogonal p functions')
        ## initialize the P matrix
        P = np.zeros((N, K))
        P[:,0] = 1.
        ## loop thru k values
        for k in range(2,K+1):
            ## determine the set associated with k
            n = self.decompose_j(k-1, Nvec)
            ## find pkhat and mu
            mu = None
            for i in range(m):
                nhat = n[:]
                if nhat[i] > 0:
                    nhat[i] -= 1
                    khat = self.compose_j(nhat, Nvec) + 1
                    if khat < k:
                        mu = i
                        break
            if mu == None: raise ValueError('Unable to find previous khat set')
            pkhat = P[:,khat-1]
            xmu = self.db.x[:,mu]
            ## calculate pk first term
            temp = xmu * pkhat
            phik = sum(n)
            pk = temp[:]
            ## loop thru summation in eq 18
            for j in range(1,k):
                ## check if value is needed
                phij = sum(self.decompose_j(j-1, Nvec))
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
            ajhat = np.dot(pj, self.db.y[:,z]) / pjdot
            ranks[i] = ajhat ** 2. * pjdot
        zm.nm.zSort(ranks, order, ascend=False, msg='Sorting the p functions by effectivenss', verbose=verbose)
        Pordered = np.zeros((N,K))
        for i,o in enumerate(order):
            Pordered[:,i] = P[:,o]
        P = Pordered[:,:]
        ###################################################################
        ##          determine how many of the orthogonal p functions to use
        if verbose: prog = zm.io.oneLineProgress(K, msg='Determining number of p functions to use')
        PSEold = None
        foundMin = False
        if sigma == None:
            yavg = sum(self.db.y[:,z]) / N
            sigma = sum([(i - yavg)**2. for i in self.db.y[:,z]]) / N
        sigma *= sigmaMultiplier
        
        for n in range(1,K+1):
            Phat = P[:,:n]
            ahat = np.matmul(np.matmul(np.linalg.inv(np.matmul(Phat.transpose(), Phat)), Phat.transpose()), self.db.y[:,z])
            yhat = np.matmul(Phat, ahat)
            MSE = np.dot(self.db.y[:,z] - yhat, self.db.y[:,z] - yhat) / N
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
        if verbose: prog = zm.io.oneLineProgress(4+nn, msg='Determining final coefficients and polynomial size')
        b = np.zeros((nn,nn))
        for k in range(1,nn+1):
            j = k - 1
            pj = P[:,j]
            w = np.ones((N,k))
            for i in range(k):
                n = self.decompose_j(order[i], Nvec)
                for ii,e in enumerate(n):
                    w[:,i] *= self.db.x[:,ii] ** e
            vals = np.matmul(np.matmul(np.linalg.inv(np.matmul(w.transpose(), w)), w.transpose()), pj)
            b[j,:k] = vals[:]
            if verbose: prog.display()
        
        A = [np.dot(P[:,i],self.db.y[:,z])/np.dot(P[:,i],P[:,i]) for i in range(nn)]
        if verbose: prog.display()
        
        c = [np.dot(A,b[:,i]) for i in range(nn)]
        js = [self.decompose_j(order[i], Nvec) for i in range(nn)]
        if verbose: prog.display()
        
        js = [js[i] for i in range(nn) if not zm.nm.isClose(c[i], 0., tol=tol)]
        c = [i for i in c if not zm.nm.isClose(i, 0., tol=tol)]
        if verbose: prog.display()
        
        nvec = [None]*m
        for i in range(m):
            nvec[i] = max([j[i] for j in js])
        JJ = 1
        for n in nvec:
            JJ *= n+1
        
        a = [0.] * JJ
        self.Jtilde[z] = 0
        for j in range(JJ):
            n = self.decompose_j(j, nvec)
            for i in range(len(c)):
                if n == js[i]:
                    a[j] = c[i]
                    self.Jtilde[z] += 1
        if verbose: prog.display()
        
        self.coef[z] = a[:]
        self.Nvec[z] = nvec[:]
        
        self.goodnessParams(mp, (z,), verbose=verbose)
    
    ########################################################################
    ########################################################################
    ########################################################################
    
    def write2Files(self):
        
        baseDir = os.path.join(os.getcwd(), self.db.name)
        
        if os.path.isdir(baseDir):
            ans = input('Data already exists. Overwrite data (y/n): ').lower()
            if ans == 'n': return
            shutil.rmtree(baseDir)
        
        os.mkdir(baseDir)
        
        for z in range(self.db.numDepVar):
            
            kw
            
            d = {
                'Independent Variable Order': self.db.namesX,
                'Nvec': self.Nvec[z],
                'coefficients': self.coef[z],
                'R2': self.R2[z],
                'RMS': self.RMS[z],
                'RMSN': self.RMSN[z],
                'Syx': self.Syx[z],
                'ybar': self.ybar[z],
                'St': self.St[z],
                'Sr': self.Sr[z],
                'settings': {k: self.kw[z][k] for k in self.kw[z] if k != 'weighting'},
                'DOF': self.Jtilde[z]
                }
            
            f = open(os.path.join(baseDir, '{}_{}.json'.format(z, self.db.namesY[z])), 'w')
            json.dump(d, f, indent=4)
            f.close()
        
        np.save(os.path.join(baseDir, 'IndependentVariables'), self.db.x)
        np.save(os.path.join(baseDir, 'DependentVariables'), self.db.y)
    

def readPolyFitsFromFiles(self, base):
    
    workingDir = os.getcwd()
    baseDir = os.path.join(workingDir, base)
    
    if not os.path.isdir(baseDir): raise ValueError()
    
    os.chdir(baseDir)
    
    kw = []
    
    
    for fn in os.listdir():
        
        if fn[-4:] == '.npy':
            
            if fn == 'DependentVariables.npy':
                x = np.load(fn)
            elif fn == 'IndependentVariables.npy':
                y = np.load(fn)
            
        elif fn[-5:] == '.json':
            
            pass
    
    
    












############################################################################
############################################################################
############################################################################

if __name__ == '__main__':
    
    k = 1000
    V = 5
    U = 10
    x = np.random.rand(k,V)
    y = np.random.rand(k,U)
    
    db = database(x, y)
    
    d1 = {
        'Nvec': [2,1,1,1,1]
        }
    
    d2 = {
        'Nvec': [1,2,1,1,1]
        }
    
    d3 = {'maxOrder': 2}
    
    kw = [d1, d2, d3, d1, d2, d1, d2, d1, d2, d1]
    
    myFits = polyFit(db, kw)
    
