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
    
    def plotSnapshot1var(self, ax, constraints, iy, f=None, avgLines=True, tol=1e-6, wireFrameColors=None, view=[30.]*2, thinning=None, **kwargsScatter):
        
        if constraints.count(None) != 2: raise ValueError('There needs to be 2 non-constraints corresponding to the two horizontal axis on the plot')
        
        I, J = -1, -1
        for i,c in enumerate(constraints):
            if I == -1 and c == None:
                I = i
            elif J == -1 and c == None:
                J = i
        
        Ncon = len(constraints)
        
        ## collect all the points that meet the constraints
        plotx, ploty, plotz = [], [], []
        if type(f) != type(None): plotf = []
        for i in range(self.numPoints):
            addPoint = True
            for j in range(Ncon):
                if constraints[j] == None: continue
                if not zm.nm.isClose(self.x[i,j], constraints[j], tol=tol):
                    addPoint = False
                    break
            if addPoint:
                plotx.append( self.x[i,I] )
                ploty.append( self.x[i,J] )
                plotz.append( self.y[i,iy] )
                if type(f) != type(None): plotf.append( f[i] )
        
        ## find unique points along the two independent variable directions
        ux, uy = [], []
        n = len(plotx)
        for i in plotx:
            flag = True
            for j in ux:
                if zm.nm.isClose(i, j, tol=tol): flag = False
            if flag: ux.append(i)
        for i in ploty:
            flag = True
            for j in uy:
                if zm.nm.isClose(i, j, tol=tol): flag = False
            if flag: uy.append(i)
        
        ## sort the unique point arrays
        zm.nm.zSort(ux, verbose=False)
        zm.nm.zSort(uy, verbose=False)
        i = len(ux)
        j = len(uy)
        
        ## thin the data if needed
        if thinning != None:
            if i > thinning:
                diff = i - thinning
                if diff % 2 == 1: diff -= 1
                # if diff == 0: diff = 2
                diff += 2
                rpX = [int(round(ind)) for ind in np.linspace(i-1,0,diff)][1:-1]
                for ind in rpX: ux.pop(ind)
            if j > thinning:
                diff = j - thinning
                if diff % 2 == 1: diff -= 1
                diff += 2
                rpY = [int(round(ind)) for ind in np.linspace(j-1,0,diff)][1:-1]
                for ind in rpY: uy.pop(ind)
            
            ## remove the thinned out points from plot arrays
            for i in range(len(plotx)-1,-1,-1):
                flag = [True, True]
                for j in range(len(ux)):
                    if zm.nm.isClose(plotx[i], ux[j], tol=tol):
                        flag[0] = False
                        break
                for j in range(len(uy)):
                    if zm.nm.isClose(ploty[i], uy[j], tol=tol):
                        flag[1] = False
                        break
                if flag[0] or flag[1]:
                    plotx.pop(i)
                    ploty.pop(i)
                    plotz.pop(i)
                    if type(f) != type(None): plotf.pop(i)
        
        n = len(plotx)
        
        ## setup meshes
        i = len(ux)
        j = len(uy)
        
        xmesh, ymesh = np.meshgrid(ux, uy)
        zmesh = np.array([[None for _ in range(i)] for _ in range(j)], dtype=float)
        if type(f) != type(None): fmesh = np.array([[None for _ in range(i)] for _ in range(j)], dtype=float)
        
        def myIndex(l, v, tol=1e-12):
            for i,j in enumerate(l):
                if zm.nm.isClose(j, v, tol=tol): return i
        
        ## fill in meshes
        for i in range(n):
            px, py, pz = plotx[i], ploty[i], plotz[i]
            row = myIndex(uy, py, tol=tol)
            col = myIndex(ux, px, tol=tol)
            # if None in (row, col): continue
            zmesh[row, col] = pz
            if type(f) != type(None): fmesh[row, col] = plotf[i]
        
        if hasattr(ax, 'get_zlim'):
            
            ## plot the wireframes
            if wireFrameColors != None:
                ax.plot_wireframe(xmesh, ymesh, zmesh, colors=wireFrameColors)
                if type(f) != type(None): ax.plot_wireframe(xmesh, ymesh, fmesh, colors=wireFrameColors)
            else:
                ax.plot_wireframe(xmesh, ymesh, zmesh, color='C0')
                if type(f) != type(None): ax.plot_wireframe(xmesh, ymesh, fmesh, color='C1')
            
            ## plot the scatter points
            ax.scatter(plotx, ploty, plotz, c='C0', **kwargsScatter)
            if type(f) != type(None): ax.scatter(plotx, ploty, plotf, c='C1', **kwargsScatter)
            ax.set_xlabel(self.namesX[I])# labelsI[I])
            ax.set_ylabel(self.namesX[J])
            ax.set_zlabel(self.namesY[iy])
            
            ## put on the avg lines
            if avgLines:
                if len(plotx) != 0:
                    avgx = sum(plotx) / len(plotx)
                    minx = min(plotx)
                    maxx = max(plotx)
                else:
                    avgx = 0.
                    minx = 0.
                    maxx = 0.
                if len(ploty) != 0:
                    avgy = sum(ploty) / len(ploty)
                    miny = min(ploty)
                    maxy = max(ploty)
                else:
                    avgy = 0.
                    miny = 0.
                    maxy = 0.
                if len(plotz) != 0:
                    avgz = sum(plotz) / len(plotz)
                    minz = min(plotz)
                    maxz = max(plotz)
                else:
                    avgz = 0.
                    minz = 0.
                    maxz = 0.
                
                ax.plot([minx,maxx], [avgy,avgy], [avgz,avgz], 'r')
                ax.plot([avgx,avgx], [miny,maxy], [avgz,avgz], 'r')
                ax.plot([avgx,avgx], [avgy,avgy], [minz,maxz], 'r')
            
            ## update the view angle
            ax.view_init(*view)
            
        else:
            
            levels = [-20, -15, -10, -5, 0, 5, 10, 15, 20]
            
            zctr = ax.contourf(xmesh, ymesh, zmesh, levels, linestyles='solid', **kwargsScatter)
            # zcbar = zm.plt.colorbar(zctr, ax=ax)
            # zcbar.ax.set_ylabel(self.namesY[iy])
            ax.clabel(zctr)
            
            ax.set_xlabel(self.namesX[I])# labelsI[I])
            ax.set_ylabel(self.namesX[J])
            ax.set_title(self.namesY[iy])
            
            if type(f) != type(None):
                fctr = ax.contour(xmesh, ymesh, fmesh, zctr.levels, linestyles='dashdot', alpha=0.8)
                # fcbar = zm.plt.colorbar(fctr, ax=ax)
                # ax.clabel(fctr)
            
            
    
    def plotSnapshot(self, fig, ax, constraints, f=None, avgLines=True, tol=1e-6, wireFrameColors=None, spa={}, view=[30.]*2, thinning=None, **kwargsScatter):
        if wireFrameColors == None: wireFrameColors = [None]*self.numDepVar
        if type(f) == type(None):
            for i in range(self.numDepVar): self.plotSnapshot1var(ax[i], constraints, i, avgLines=avgLines, tol=tol, wireFrameColors=wireFrameColors[i], view=view, thinning=thinning, **kwargsScatter)
        else:
            F = np.asarray(f)
            for i in range(self.numDepVar): self.plotSnapshot1var(ax[i], constraints, i, f=F[:,i], avgLines=avgLines, tol=tol, wireFrameColors=wireFrameColors[i], view=view, thinning=thinning, **kwargsScatter)
        numConstVar = self.numIndVar - 2
        ii = [i for i,j in enumerate(constraints) if j != None]
        C = [self.namesX[i] for i in ii]
        vals = [constraints[i] for i in ii]
        fig.suptitle(('  {} = {}'*numConstVar).format(*[j for i in zip(C, vals) for j in i]))
        if spa != {}:
            fig.subplots_adjust(**spa)
        else:
            fig.tight_layout()
    
    def viewData(self, fig, ax, f=None, wireFrameColors=None, spa={}, tol=1e-6, **kwargsScatter):
        
        if wireFrameColors == None: wireFrameColors = [None]*self.numDepVar
        
        prec = str(int(abs(np.floor(np.log10(tol)))))
        
        numConstVar = self.numIndVar - 2
        
        cont = True
        while cont:
            
            C = []
            if numConstVar > 0: zm.io.text(*self.namesX, title='Choose {} variables to hold constant'.format(numConstVar))
            while len(C) != numConstVar:
                incorrect = True
                while incorrect:
                    c = input('Enter variable:  ')
                    if c in self.namesX and c not in C:
                        incorrect = False
                    else:
                        print('invalid entry, try again')
                C.append(c)
            print()
            
            ii = [self.namesX.index(c) for c in C]
            Consts = [None]*self.numIndVar
            vals = [None]*numConstVar
            
            for i in range(numConstVar):
                
                ## find unique points along the independent variable direction
                u = []
                for I in self.x[:,ii[i]]:
                    flag = True
                    for j in u:
                        if zm.nm.isClose(I, j, tol=tol): flag = False
                    if flag: u.append(I)
                ## sort the unique point arrays
                zm.nm.zSort(u, verbose=False)
                
                print('Unique {} values in database:'.format(C[i]))
                print((('  {:.'+prec+'f}')*len(u)).format(*u))
                vals[i] = float(input('Choose a value for {}: '.format(C[i])))
                Consts[ii[i]] = vals[i]
                print()
            
            if hasattr(ax[0], 'get_zlim'):
                
                ele = float(input('Enter elevation view angle for the plot(s) in degrees: '))
                rot = float(input('Enter rotaion view angle for the plot(s) in degrees:   '))
                print()
                
                thin = int(input('Enter desired number of points along each dimension to plot, 0 for all points: '))
                
            else:
                ele = rot = thin = 0
            if thin == 0: thin = None
            # axTemp = fig.axes
            # print(len(ax), len(axTemp))
            # for i in range(len(ax),len(axTemp)): axTemp[i].remove()
            
            # print(len(fig.axes))
            
            for i in ax: i.cla()
            
            self.plotSnapshot(fig, ax, Consts, f=f, wireFrameColors=wireFrameColors, spa=spa, view=[ele, rot], thinning=thin, **kwargsScatter)
            
            # fig.suptitle(('  {} = {}'*numConstVar).format(*[j for i in zip(C, vals) for j in i]))#C[0], vals[0], C[1], vals[1]))
            
            if input('Plot again (y/n)? ').lower() == 'n': cont = False
            print()
    

class polyFit():
    
    def __init__(self, *args, mpFits=1, verbose=True, c=0):
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
            independent variables and can be performed with full control of
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
                    list of tuples containing two integers. The two integers
                    represent the independent variables that the even cross
                    symmetry is applied, which is an even symmetry that acts
                    along a diagonal among the two independent variables 
                    listed. The even cross symmetry forces all polynomial
                    terms between the two independent variables that have
                    different parity on the exponents to zero, i.e. x1^4 *
                    x2^5 have exponents that are odd and even, so the even
                    cross symmetry would force the poly coefficient related
                    to this term to 0. Defaults to an empty list
                
                crossSymOdd : list, optional
                    list of tuples containing two integers. The two integers
                    represent the independent variables that the odd cross
                    symmetry is applied, which is an odd symmetry that acts
                    along a diagonal among the two independent variables 
                    listed. The odd cross symmetry forces all polynomial
                    terms between the two independent variables that have
                    the same parity on the exponents to zero, i.e. x1^2*x2^4
                    both exponents are even so the odd cross symmetry would
                    force the poly coefficient related to this term to 0.
                    Defaults to an empty list
                
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
        
        ## check if number of inputs is one
        if len(args) == 1:
            self.readPolyFitsFromFiles(args[0], verbose=verbose)
            self.c = c
        else:
            db, kw = args
            
            ## copy in database
            self.db = db
            
            self.c = c
            
            ## create array for polynomial values corresponding to the db values
            self.f = np.zeros((self.db.numPoints, self.db.numDepVar)) * np.nan
            
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
            self.Jtilde  = np.zeros( self.db.numDepVar, dtype=int )
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
                    if self.kw[unusedManFits[i]] == self.kw[duplicates[0]] and not self.kw[duplicates[0]].get('percent', False) and not 'weighting' in self.kw[duplicates[0]]:
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
                    if verbose: zm.io.text('Performing manual fit(s) for:', *['{}'.format(self.db.namesY[i]) for i in dupFit], c=self.c)
                    self.manFit(dupFit)
                ## peform the auto fits
                for i in autoFits:
                    if verbose: zm.io.text('Performing auto fit for {}'.format(self.db.namesY[i]), c=self.c)
                    self.autoFit(i)
            else: ## perfoming the fits simultanuously with multiprocessing
                ## disable the multiprocessing option for all fits and verbosity
                for i in range(self.db.numDepVar):
                    self.kw[i]['mp'] = 1
                    self.kw[i]['verbose'] = False
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
                if verbose: prog = zm.io.Progress(len(it), title='Performing fits for {}'.format(self.db.name), c=self.c)
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
    
    def computeWeighting(self, args):
        kk, z, weighting, percent = args
        w = 1.
        if callable(weighting): w *= weighting(self.db, z, kk) ** 2.
        if percent: w /= self.db.y[kk, z] **2.
        return kk, w
    
    def __computeCHU__(self, n, cpus):
        chu = n // cpus // 20
        # chu = 10
        if chu > 8000: return 8000
        if chu < 1: return 1
        return chu
    
    def createA(self, args):
        r, c, z, w, p = args
        a = 0.
        nr = self.decompose_j(r, self.Nvec[z])
        nc = self.decompose_j(c, self.Nvec[z])
        n = [nr[i] + nc[i] for i in range(self.db.numIndVar)]
        for kk in range(self.db.numPoints):
            t = self.computeWeighting((kk, z, w, p))[-1]
            for v in range(self.db.numIndVar):
                for _ in range(n[v]): t *= self.db.x[kk,v]
            a += t
        return r, c, a
    
    def createB(self, args):
        r, iy, w, p = args
        b = np.zeros(len(iy))
        nr = self.decompose_j(r, self.Nvec[iy[0]])
        for ib,z in enumerate(iy):
            for kk in range(self.db.numPoints):
                t = self.computeWeighting((kk, z, w, p))[-1] * self.db.y[kk,z]
                for v in range(self.db.numIndVar):
                    for _ in range(nr[v]): t *= self.db.x[kk,v]
                b[ib] += t
        return r, b
    
    def manFit(self, iy):
        
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
        saveMemory      = self.kw[z].get('saveMemory', False)
        
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
            ## loop through crossSymOdd constraints
            for val in crossSymOdd:
                if flag: break
                ## check if the n values from both variables given in val are even, then trip flag
                if n[val[0]]%2 == 0 and n[val[1]]%2 == 0:
                    flag = True
                ## check if the n values from both variables given in val are odd, then trip flap
                if n[val[0]]%2 == 1 and n[val[1]]%2 == 1:
                    flag = True
            if flag: continue
            ## loop through crossSymEven constraints
            for val in crossSymEven:
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
        
        ## compute A and b matrices
        ########################################################################
        ## attempt to use the basis functions matrix
        try:
            if saveMemory: raise MemoryError
            X = np.zeros((k,lenActive))
            if callable(weighting) or percent: Xt = X.T.copy()
            A = np.zeros((lenActive,lenActive))
            b = np.zeros((lenActive, len(iy)))
            
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
                    for vals in pool.imap_unordered(self.createX, it, chunksize=self.__computeCHU__(lenActive*k, cpus)):
                        kk, jj, val = vals
                        X[kk, jj] = val
                        if verbose: prog.display()
                del it
            
            if callable(weighting) or percent:
                Xt = X.T.copy()
                if verbose: prog = zm.io.oneLineProgress(k, msg='Setting up weighting factors')
                if mp == 1:
                    for kk in range(k):
                        Xt[:,kk] *= self.computeWeighting((kk, z, weighting, percent))[-1]
                        if verbose: prog.display()
                else:
                    it = [(kk, z, weighting, percent) for kk in range(k)]
                    with Pool(cpus) as pool:
                        for kk, w in pool.imap_unordered(self.computeWeighting, it, chunksize=self.__computeCHU__(k, cpus)):
                            Xt[:,kk] *= w
                            if verbose: prog.display()
                    del it
                
                if verbose: zm.io.oneLineText('Computing the A matrix and b vector', c=self.c)
                
                A = Xt.dot(X)
                b = Xt.dot(self.db.y)
                
                del Xt, X
                
            else:
                
                if verbose: zm.io.oneLineText('Computing the A matrix and b vector', c=self.c)
                
                A = X.T.dot(X)
                b = X.T.dot(self.db.y)
                
                del X
            
        except MemoryError: ## basis function matrix method uses too much memory, must compute A and b manually
        # else:
            
            if verbose:
                zm.io.oneLineText('Insufficient memory to use Basis Function Matrix', c=self.c)
                prog = zm.io.oneLineProgress(lenActive**2+lenActive, msg='Computing the A matrix and b vector')
            
            A = np.zeros((lenActive,lenActive))
            b = np.zeros((lenActive,len(iy)))
            
            if mp == 1:
                
                for row in range(lenActive):
                    for col in range(lenActive):
                        A[row,col] = self.createA((row,col,z,weighting,percent))[-1]
                        if verbose: prog.display()
                    b[row,:] = self.createB(row,iy,weighting,percent)[-1]
                    if verbose: prog.display()
                
            else:
                if mp == 0:
                    cpus = cpu_count()
                else:
                    cpus = mp
                
                itA = [None]*(lenActive**2)
                itB = [None]*lenActive
                i = -1
                for row in range(lenActive):
                    for col in range(lenActive):
                        i += 1
                        itA[i] = (row, col, z, weighting, percent)
                    itB[row] = (row, iy, weighting, percent)
                
                with Pool(cpus) as pool:
                    for row, col, val in pool.imap_unordered(self.createA, itA, chunksize=self.__computeCHU__(lenActive**2, cpus)):
                        A[row,col] = val
                        if verbose: prog.display()
                    for row, val in pool.imap_unordered(self.createB, itB, chunksize=self.__computeCHU__(lenActive, cpus)):
                        b[row,:] = val
                        if verbose: prog.display()
                del itA, itB
        
        
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
        if verbose: zm.io.oneLineText('solving the Aa=b equation', c=self.c)
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
            ## compute weighting terms
            W = self.kw[z].get('weighting', None)
            P = self.kw[z].get('percent', False)
            w = np.copy([self.computeWeighting((kk, z, W, P))[-1] for kk in range(self.db.numPoints) if self.db.y[kk,z] != None])
            ## calculate the SSt value
            self.St[z] = float(sum( ((ynew - self.ybar[z])*w) ** 2. ))
            
            ## loop through the datapoints
            if verbose: prog = zm.io.oneLineProgress(k, msg='Evaluating Fit Parameters for {}'.format(self.db.namesY[z]))
            
            if mp == 1:
                for i in range(self.db.numPoints):
                    if self.db.y[i,z] == None:
                        continue
                    ## calculate the f value from the polynomial function
                    self.f[i,z] = self.evaluate(z, self.db.x[i,:])
                    if verbose: prog.display()
            else:
                it = []
                for i in range(self.db.numPoints):
                    if self.db.y[i,z] == None: continue
                    it.append((i,z))
                
                if mp == 0:
                    cpus = cpu_count()
                else:
                    cpus = mp
                
                with Pool(cpus) as pool:
                    chu = k // cpus // 20
                    if chu > 8000: chu = 8000
                    if chu < 1: chu = 1
                    for i, val in pool.imap_unordered(self.evalMP, it, chunksize=chu):
                        self.f[i,z] = val
                        if verbose: prog.display()
            
            fnew = np.copy([temp for temp in self.f[:,z] if temp != np.nan])
            
            
            # calculate the SSr term
            self.Sr[z] = float(sum( ((ynew - fnew)*w) ** 2. ))
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
            self.RMSN[z] = np.sqrt(np.mean(((ynew - fnew)*w/avg) ** 2.))
    
    def evalMP(self, args):
        return args[0], self.evaluate(args[1], self.db.x[args[0],:])
    
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
                for _ in range(n[v]): prod *= x[v]
                # prod *= x[v] ** n[v]
            # add onto the summation the proper term
            f += self.coef[z][j] * prod
        # return the finalized summation value
        return f
    
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
        return zm.io.text(msg, p2s=False, title=self.db.namesY[z], c=self.c)
    
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
        zm.nm.zSort(ranks, order, ascend=False, msg='Sorting the p functions by effectivenss', verbose=verbose, c=self.c)
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
                'DOF': int(self.Jtilde[z])
                }
            
            f = open(os.path.join(baseDir, '{}_{}.json'.format(z, self.db.namesY[z])), 'w')
            json.dump(d, f, indent=4)
            f.close()
        
        np.save(os.path.join(baseDir, 'IndependentVariables'), self.db.x)
        np.save(os.path.join(baseDir, 'DependentVariables'), self.db.y)
        np.save(os.path.join(baseDir, 'PolyFuncValues'), self.f)
    
    
    def readPolyFitsFromFiles(self, base, verbose=True):
        
        workingDir = os.getcwd()
        baseDir = os.path.join(workingDir, base)
        
        if not os.path.isdir(baseDir): raise ValueError()
        
        os.chdir(baseDir)
        
        order = []
        self.Nvec = []
        self.coef = []
        self.R2 = []
        self.RMS = []
        self.RMSN = []
        self.Syx = []
        self.ybar = []
        self.St = []
        self.Sr = []
        self.kw = []
        self.Jtilde = []
        namesY = []
        self.numCoef = []
        
        if verbose: prog = zm.io.oneLineProgress(len(os.listdir()), msg='Reading in files')
        
        for fn in os.listdir():
            
            if fn[-4:] == '.npy':
                
                if fn == 'IndependentVariables.npy':
                    x = np.load(fn)
                elif fn == 'DependentVariables.npy':
                    y = np.load(fn)
                elif fn == 'PolyFuncValues.npy':
                    self.f = np.load(fn)
                
            elif fn[-5:] == '.json':
                
                for I,c in enumerate(fn):
                    if c == '_': break
                
                i = int(fn[:I])
                
                order.append(i)
                
                f = open(fn, 'r')
                data = json.load(f)
                f.close()
                
                self.Nvec.append( data['Nvec'] )
                self.coef.append( data['coefficients'] )
                self.R2.append( data['R2'] )
                self.RMS.append( data['RMS'] )
                self.RMSN.append( data['RMSN'] )
                self.Syx.append( data['Syx'] )
                self.ybar.append( data['ybar'] )
                self.St.append( data['St'] )
                self.Sr.append( data['Sr'] )
                self.kw.append( data['settings'] )
                self.Jtilde.append( data['DOF'] )
                self.numCoef.append( self.calcNumCoef(data['Nvec']) )
                
                namesY.append( fn[I+1:-5] )
                
                if i == 0: namesX = data['Independent Variable Order']
            
            if verbose: prog.display()
        
        os.chdir(workingDir)
        
        zm.nm.zSort(order, self.Nvec, self.coef, self.R2, self.RMS, self.RMSN, self.Syx, self.ybar, self.St, self.Sr, self.kw, self.Jtilde, namesY, self.numCoef, verbose=verbose)
        
        self.db = database(x, y, name=base, namesX=namesX, namesY=namesY)
        
        self.auto = [False if 'Nvec' in self.kw[i] else True for i in range(self.db.numDepVar)]
    
    
    
    












############################################################################
############################################################################
############################################################################

if __name__ == '__main__':
    
    k = 10000
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
    
    d3 = {'maxOrder': 2, 'sigmaMultiplier':1e-5}
    
    # kw = [d1, d2, d3, d1, d2, d1, d2, d1, d2, d1]
    kw = [d3, d3, d3, d3, d3, d3, d3, d3, d3, d3]
    
    myFits = polyFit(db, kw, mpFits=0)
    
    for c in myFits.coef:
        print()
        print(c)
