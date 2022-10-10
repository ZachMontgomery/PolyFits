<<<<<<< HEAD
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
    https://doi.org/10.2514/3.46712

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
=======
>>>>>>> da763d0cd7e0f1fddb43f2eb2bf81427119cd3c7
import numpy as np
import ZachsModules as zm
from multiprocessing import Pool, cpu_count
import os
import shutil
import json
from matplotlib import cm

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
    
    def plotSnapshot1var(self, ax, constraints, iy, f=None, avgLines=True, tol=1e-6, wireFrameColors=None, view=[30.]*2, thinning=None, makeScatter=True, numClusters=None, **kwargsScatter):
        
        I, x, y, F = self.constrainData(constraints,tol,F=f)
        
        z = np.array(y[:,iy])
        y = np.array(x[:,1])
        x = np.array(x[:,0])
        
        xmesh, ymesh, zmesh = self.clusterMesh(x,y,z,numClusters=numClusters)
        fmesh = self.clusterMesh(x,y,F[:,0],numClusters=numClusters)[-1]
        
        if hasattr(ax, 'get_zlim'):
            
            ## plot the wireframes
            if wireFrameColors != None:
                MESH = ax.plot_wireframe(xmesh, ymesh, zmesh, colors=wireFrameColors)
                if type(f) != type(None): ax.plot_wireframe(xmesh, ymesh, fmesh, colors=wireFrameColors)
            else:
                if type(f) == type(None):
                    norm = zm.plt.Normalize(*ax.get_zlim3d())
                    colors = cm.viridis(norm(zmesh))
                    MESH = ax.plot_surface(xmesh, ymesh, zmesh, facecolors=colors, shade=False, linewidth=0.5)
                    MESH.set_facecolor((0,0,0,0))
                    ax.contour(xmesh, ymesh, zmesh, zdir='z', offset=ax.get_zlim3d()[0], cmap=cm.viridis, vmin=ax.get_zlim3d()[0], vmax=ax.get_zlim3d()[1])
                if type(f) != type(None):
                    MESH = ax.plot_wireframe(xmesh, ymesh, zmesh, color='C0')
                    ax.plot_wireframe(xmesh, ymesh, fmesh, color='C1')
            
            ## plot the scatter points
            if makeScatter:
                if type(f) == type(None):
                    ax.scatter(x, y, z, c='r', **kwargsScatter)
                if type(f) != type(None):
                    ax.scatter(x, y, z, c='C0', **kwargsScatter)
                    ax.scatter(x, y, F, c='C1', **kwargsScatter)
            ax.set_xlabel(self.namesX[I[0]])# labelsI[I])
            ax.set_ylabel(self.namesX[I[1]])
            ax.set_zlabel(self.namesY[iy])
            
            ## put on the avg lines
            if avgLines:
                if len(x) != 0:
                    avgx = sum(x) / len(x)
                    minx = min(x)
                    maxx = max(x)
                else:
                    avgx = 0.
                    minx = 0.
                    maxx = 0.
                if len(y) != 0:
                    avgy = sum(y) / len(y)
                    miny = min(y)
                    maxy = max(y)
                else:
                    avgy = 0.
                    miny = 0.
                    maxy = 0.
                if len(z) != 0:
                    avgz = sum(z) / len(z)
                    minz = min(z)
                    maxz = max(z)
                else:
                    avgz = 0.
                    minz = 0.
                    maxz = 0.
                
                ax.plot([minx,maxx], [avgy,avgy], [avgz,avgz], 'r')
                ax.plot([avgx,avgx], [miny,maxy], [avgz,avgz], 'r')
                ax.plot([avgx,avgx], [avgy,avgy], [minz,maxz], 'r')
            
            ## update the view angle
            ax.view_init(*view)
            
            return MESH
        else:
            
            # levels = [-20, -15, -10, -5, 0, 5, 10, 15, 20]
            
            zctr = ax.contourf(xmesh, ymesh, zmesh, linestyles='solid', **kwargsScatter)
            # zcbar = zm.plt.colorbar(zctr, ax=ax)
            # zcbar.ax.set_ylabel(self.namesY[iy])
            ax.clabel(zctr)
            
            ax.set_xlabel(self.namesX[I[0]])# labelsI[I])
            ax.set_ylabel(self.namesX[I[1]])
            ax.set_title(self.namesY[iy])
            
            if type(f) != type(None):
                fctr = ax.contour(xmesh, ymesh, fmesh, zctr.levels, linestyles='dashdot', alpha=0.8)
                # fcbar = zm.plt.colorbar(fctr, ax=ax)
                # ax.clabel(fctr)
            
            return zctr
    
    def plotSnapshot1varOld(self, ax, constraints, iy, f=None, avgLines=True, tol=1e-6, wireFrameColors=None, view=[30.]*2, thinning=None, makeScatter=True, **kwargsScatter):
        
        if constraints.count(None) != 2: raise ValueError('There needs to be 2 non-constraints corresponding to the two horizontal axis on the plot')
        
        if not zm.misc.isIterable(tol): tol = [tol]*self.numIndVar
        
        I, J = None, None
        for i,c in enumerate(constraints):
            if I == None and c == None:
                I = i
            elif J == None and c == None:
                J = i
                break
        
        Ncon = len(constraints)
        
        ## collect all the points that meet the constraints
        plotx, ploty, plotz = [], [], []
        if type(f) != type(None): plotf = []
        for i in range(self.numPoints):
            addPoint = True
            for j in range(Ncon):
                if constraints[j] == None: continue
                if not zm.nm.isClose(self.x[i,j], constraints[j], tol=tol[j]):
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
                if zm.nm.isClose(i, j, tol=tol[I]): flag = False
            if flag: ux.append(i)
        for i in ploty:
            flag = True
            for j in uy:
                if zm.nm.isClose(i, j, tol=tol[J]): flag = False
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
                    if zm.nm.isClose(plotx[i], ux[j], tol=tol[I]):
                        flag[0] = False
                        break
                for j in range(len(uy)):
                    if zm.nm.isClose(ploty[i], uy[j], tol=tol[J]):
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
            row = myIndex(uy, ploty[i], tol=tol[J])
            col = myIndex(ux, plotx[i], tol=tol[I])
            zmesh[row, col] = plotz[i]
            if type(f) != type(None): fmesh[row, col] = plotf[i]
        
        if hasattr(ax, 'get_zlim'):
            
            ## plot the wireframes
            if wireFrameColors != None:
                MESH = ax.plot_wireframe(xmesh, ymesh, zmesh, colors=wireFrameColors)
                if type(f) != type(None): ax.plot_wireframe(xmesh, ymesh, fmesh, colors=wireFrameColors)
            else:
                if type(f) == type(None):
                    norm = zm.plt.Normalize(*ax.get_zlim3d())
                    colors = cm.viridis(norm(zmesh))
                    MESH = ax.plot_surface(xmesh, ymesh, zmesh, facecolors=colors, shade=False, linewidth=0.5)
                    MESH.set_facecolor((0,0,0,0))
                    ax.contour(xmesh, ymesh, zmesh, zdir='z', offset=ax.get_zlim3d()[0], cmap=cm.viridis, vmin=ax.get_zlim3d()[0], vmax=ax.get_zlim3d()[1])
                if type(f) != type(None):
                    MESH = ax.plot_wireframe(xmesh, ymesh, zmesh, color='C0')
                    ax.plot_wireframe(xmesh, ymesh, fmesh, color='C1')
            
            ## plot the scatter points
            if makeScatter:
                if type(f) == type(None):
                    ax.scatter(plotx, ploty, plotz, c='r', **kwargsScatter)
                if type(f) != type(None):
                    ax.scatter(plotx, ploty, plotz, c='C0', **kwargsScatter)
                    ax.scatter(plotx, ploty, plotf, c='C1', **kwargsScatter)
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
            
            return MESH
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
            
            return zctr
    
    def constrainData(self,constraints,tol,F=None):
        if not zm.misc.isIterable(tol): tol = [tol]*self.numIndVar
        if constraints.count(None) != 2 or len(constraints) != self.numIndVar: raise ValueError('There needs to be 2 non-constraints corresponding to the two horizontal axis on the plot and a constraint for each remaining independent variable.')
        if F != None:
            f = np.array(F)
        else:
            f = np.zeros(self.y.shape)
        ## get indices of the two independent variables
        I = []
        K = []
        for i,c in enumerate(constraints):
            if c == None:
                I.append(i)
            else:
                K.append(i)
        ## initializations
        x = np.array(self.x)
        y = np.array(self.y)
        x = list(x)
        x = [list(i) for i in x]
        y = list(y)
        y = [list(i) for i in y]
        f = list(f)
        f = [list(i) for i in f]
        ## remove points that don't meet the constraints within the tolerances
        for i in range(len(x)-1,-1,-1):
            meetsCon = True
            for v in K:
                if not zm.nm.isClose(x[i][v], constraints[v], tol=tol[v]):
                    meetsCon = False
                    break
            if not meetsCon:
                x.pop(i)
                y.pop(i)
                f.pop(i)
        ## return results
        x = np.array(x)
        y = np.array(y)
        f = np.array(f)
        return I, x[:,I], y, f
    
    @staticmethod
    def closestIndex(ar,val):
        d = [abs(val - a) for a in ar]
        return d.index(min(d))
    
    @staticmethod
    def maxDiscrepency(mesh,ar):
        D = []
        for m in mesh:
            d = [abs(m-a) for a in ar]
            D += [min(d)]
        return max(D) / (max(ar) - min(ar)) * 100
    
    @staticmethod
    def clusterData2numLvls(D,n):
        
        def printgr(group):
            print()
            print([np.mean(gr) for gr in group])
            print([np.std(gr) for gr in group])
            print([len(gr) for gr in group])
            print()
        
        def regroup(mesh,data):
            group = [[] for _ in range(len(mesh))]
            for d in data:
                i = database.closestIndex(mesh,d)
                group[i].append(d)
            return group
        
        d = np.array(D)
        if n > len(d): raise ValueError('Ummm, you do not have enough datapoints, {}, for that many levels, {}.'.format(len(d),n))
        if n < 1: raise ValueError('Why would you do this? Just solve for the average if that is what you want.')
        
        ## first grouping
        mesh = [min(d),max(d)]
        
        group = regroup(mesh,d)
        
        mesh = [np.mean(gr) for gr in group]
        group = regroup(mesh,d)
        
        for _ in range(n-2):
            ## update xmesh
            std = [np.std(gr) for gr in group]
            i = std.index(max(std))
            m,M = min(group[i]),max(group[i])
            mesh.pop(i)
            mesh.append(m)
            mesh.append(M)
            zm.nm.zSort(mesh, verbose=False)
            
            group = regroup(mesh,d)
            
            mesh = [np.mean(gr) for gr in group]
            group = regroup(mesh, d)
        return mesh
    
    @staticmethod
    def createZmesh(ux, uy, X, Y, Z):
        nx, ny = len(ux), len(uy)
        
        zmesh = np.zeros((ny,nx)) * np.nan
        
        dup = {}
        
        for i in range(len(Z)):
            x = database.closestIndex(ux,X[i])
            y = database.closestIndex(uy,Y[i])
            if np.isnan(zmesh[y,x]):
                zmesh[y,x] = Z[i]
            else:
                if (y,x) in dup:
                    dup[(y,x)].append(Z[i])
                else:
                    dup[(y,x)] = [zmesh[y,x],Z[i]]
                zmesh[y,x] = np.mean(dup[(y,x)])
        
        return zmesh
    
    @staticmethod
    def clusterMesh(X,Y,Z,numClusters=None):
        '''
        Takes 1D arrays for X, Y, and Z
        returns the finest mesh grid arrays that will NOT have any gaps
        '''
        x = np.array(X)
        y = np.array(Y)
        z = np.array(Z)
        n = len(x)
        
        if  len(x.shape) > 1 | len(y.shape) > 1 | len(z.shape) > 1 | n != len(y) | n != len(z): raise ValueError()
        
        if numClusters != None:
            ux = database.clusterData2numLvls(x,numClusters[0])
            uy = database.clusterData2numLvls(y,numClusters[1])
            zmesh = database.createZmesh(ux,uy,x,y,z)
            xmesh, ymesh = np.meshgrid(ux, uy)
            return xmesh, ymesh, zmesh
        
        nx = ny = 2
        
        ux = database.clusterData2numLvls(x,nx)
        uy = database.clusterData2numLvls(y,ny)
        zmesh = database.createZmesh(ux, uy, x, y, z)
        
        if np.isnan(zmesh).sum() > 0: raise ValueError('Cannot even mesh 2x2 grid')
        
        while np.isnan(zmesh).sum() == 0:
            
            dx = database.maxDiscrepency(ux,x)
            dy = database.maxDiscrepency(uy,y)
            if dx > dy:
                incX = True
            else:
                incX = False
            if incX:
                nx += 1
            else:
                ny += 1
            ux = database.clusterData2numLvls(x,nx)
            uy = database.clusterData2numLvls(y,ny)
            zmesh = database.createZmesh(ux,uy,x,y,z)
            # print()
            # print(nx,ny)
            # print(ux)
            # print(uy)
            # print(np.isnan(zmesh).sum())
            # input()
        
        if incX:
            nx -= 1
        else:
            ny -= 1
        ux = database.clusterData2numLvls(x,nx)
        uy = database.clusterData2numLvls(y,ny)
        zmesh = database.createZmesh(ux,uy,x,y,z)
        xmesh, ymesh = np.meshgrid(ux, uy)
        
        return xmesh, ymesh, zmesh
    
    
    
    
    def plotSnapshot(self, fig, ax, iy, constraints, f=None, avgLines=True, tol=1e-6, wireFrameColors=None, spa={}, view=[30.]*2, thinning=None, makeScatter=True, numClusters=None, **kwargsScatter):
        if not zm.misc.isIterable(tol): tol = [tol]*self.numIndVar
        if wireFrameColors == None: wireFrameColors = [None]*self.numDepVar
        MESHES = [None]*self.numDepVar
        if type(f) == type(None):
            for i,j in enumerate(iy):
                MESHES[i] = self.plotSnapshot1var(ax[i], constraints, j, avgLines=avgLines, tol=tol, wireFrameColors=wireFrameColors[i], view=view, thinning=thinning, makeScatter=makeScatter, numClusters=numClusters, **kwargsScatter)
        else:
            F = np.asarray(f)
            for i,j in enumerate(iy):
                MESHES[i] = self.plotSnapshot1var(ax[i], constraints, j, f=F[:,j], avgLines=avgLines, tol=tol, wireFrameColors=wireFrameColors[i], view=view, thinning=thinning, makeScatter=makeScatter, numClusters=numClusters, **kwargsScatter)
        numConstVar = self.numIndVar - 2
        ii = [i for i,j in enumerate(constraints) if j != None]
        C = [self.namesX[i] for i in ii]
        vals = [constraints[i] for i in ii]
        tols = [tol[i] for i in ii]
        fig.suptitle((r'  {} = {}$\pm${}'*numConstVar).format(*[j for i in zip(C, vals, tols) for j in i]))
        if spa != {}:
            fig.subplots_adjust(**spa)
        else:
            fig.tight_layout()
        # fig.canvas.draw_idle()
        return MESHES
    
    def viewData(self, fig, ax, iy, f=None, wireFrameColors=None, spa={}, tol=1e-6, zlim=(), avgLines=True, makeScatter=True, numClusters=None, **kwargsScatter):
        
        if wireFrameColors == None: wireFrameColors = [None]*self.numDepVar
        
        if not zm.misc.isIterable(tol): tol = [tol]*self.numIndVar
        prec = [str(int(abs(np.floor(np.log10(i))))) for i in tol]
        
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
                    if I == None or np.isnan(I): continue
                    flag = True
                    for j in u:
                        if zm.nm.isClose(I, j, tol=tol[ii[i]]): flag = False
                    if flag: u.append(I)
                ## sort the unique point arrays
                zm.nm.zSort(u, verbose=False)
                
                print('Unique {} values in database:'.format(C[i]))
                print((('  {:.'+prec[ii[i]]+'f}')*len(u)).format(*u))
                vals[i] = float(input('Choose a value for {}: '.format(C[i])))
                Consts[ii[i]] = vals[i]
                print()
            
            if hasattr(ax[0], 'get_zlim'):
                
                ele = float(zm.io.timedInput('Enter elevation view angle for the plot(s) in degrees',15.0,timeout=1))
                rot = float(zm.io.timedInput('Enter rotaion view angle for the plot(s) in degrees  ',35.0,timeout=1))
                print()
                
                thin = int(zm.io.timedInput('Enter desired number of points along each dimension to plot, 0 for all points',0,timeout=1))
                
            else:
                ele = rot = thin = 0
            if thin == 0: thin = None
            # axTemp = fig.axes
            # print(len(ax), len(axTemp))
            # for i in range(len(ax),len(axTemp)): axTemp[i].remove()
            
            # print(len(fig.axes))
            
            for i in ax: i.cla()
            if len(zlim) == len(ax):
                for j,i in enumerate(ax):
                    if hasattr(i, 'get_zlim'):
                        i.set_zlim3d(zlim[j])
            
            meshes = self.plotSnapshot(fig, ax, iy, Consts, f=f, wireFrameColors=wireFrameColors, spa=spa, view=[ele, rot], thinning=thin, tol=tol, avgLines=avgLines, makeScatter=makeScatter, numClusters=numClusters, **kwargsScatter)
            
            # fig.suptitle(('  {} = {}'*numConstVar).format(*[j for i in zip(C, vals) for j in i]))#C[0], vals[0], C[1], vals[1]))
            
            if input('Plot again (y/n)? ').lower() == 'n': cont = False
            print()
        return meshes
    

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
            self.c = c
            self.readPolyFitsFromFiles(args[0], verbose=verbose)
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
                if verbose: prog = zm.io.oneLineProgress(len(it), msg='Performing fits for {}'.format(self.db.name), c=self.c)
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
        interaction     = self.kw[z].get('interaction', [])
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
        
        if interaction == True:
            interaction = []
        elif interaction == False:
            interaction = [[i,j] for i in range(self.db.numIndVar-1) for j in range(i,self.db.numIndVar) if i != j]
        
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
            ## check if j is an allowed interaction term and if not then continue on to the next j
            
            ## set interaction flag to true
            interactionFlag = True
            ## loop thru unallowed interactions
            for vals in interaction:
                if not 0 in [n[val] for val in vals]:
                    interactionFlag = False
                    break
            if not interactionFlag: continue
            
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
                        if n[count]%2 == 0 and n[count] != 0:
                            flag = True
            ## if the flag has been tripped, skip to the next j value
            if flag: continue
            ## loop through crossSymOdd constraints
            for val in crossSymOdd:
                if flag: break
                ## check if the n values from both variables given in val are even, then trip flag
                if n[val[0]]%2 == 0 and n[val[1]]%2 == 0:
                    if n[val[0]] == 0 and n[val[1]] == 0: continue  ## allow offset term
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
            if verbose: prog = zm.io.oneLineProgress(k*lenActive, msg='PolyFit Setup: Computing the Basis Functions', c = self.c)
            
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
                if verbose: prog = zm.io.oneLineProgress(k, msg='Setting up weighting factors', c = self.c)
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
                prog = zm.io.oneLineProgress(lenActive**2+lenActive, msg='Computing the A matrix and b vector', c = self.c)
            
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
            if verbose: prog = zm.io.oneLineProgress(k, msg='Evaluating Fit Parameters for {}'.format(self.db.namesY[z]), c = self.c)
            
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
        mp              = self.kw[z].get('mp', 0)
        
        m = self.db.numIndVar
        ## max range of polynomials to try
        Nvec = tuple([maxOrder]*m)
        ## number of datapoints
        N = self.db.numPoints
        ## number of p functions
        K = self.kCompose(Nvec)
        ###################################################################
        ##           determine the orthogonal p functions
        if verbose: prog = zm.io.oneLineProgress(K-1, msg='Determining the orthogonal p functions', c = self.c)
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
        if verbose: prog = zm.io.oneLineProgress(K, msg='Determining number of p functions to use', c = self.c)
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
        if verbose: prog = zm.io.oneLineProgress(4+nn, msg='Determining final coefficients and polynomial size', c = self.c)
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
        self.numCoef[z] = self.calcNumCoef(nvec)
        
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
        
        if verbose: prog = zm.io.oneLineProgress(len(os.listdir()), msg='Reading in files', c = self.c)
        
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
    
    
    def writeHardCodedEqs(self, **kws):
        
        namesX = kws.get('namesX', self.db.namesX)
        namesY = kws.get('namesY', self.db.namesY)
        latex  = kws.get('latex' , False)
        
        
        tab = ' '*4
        nl = '\n'
        
        s = ''
        
        for i in range(self.db.numDepVar):
            
            s += nl*2 + 'def evaluate_' + namesY[i] + '('+('{}, '*self.db.numIndVar).format(*[v+'1' for v in namesX])+'):' + nl
            
            
            # s += tab + ('{} = '*self.db.numIndVar).format(*[v+'0' for v in namesX]) + '1.0' + nl
            
            for v in range(self.db.numIndVar):
                if self.Nvec[i][v] > 1:
                    s += tab + nl
                    for j in range(2, self.Nvec[i][v]+1):
                        s += tab + '{0}{1} = {0}{2} * {0}1'.format(namesX[v], j, j-1)+nl
            
            s += tab + nl
            
            s += tab + 'return (' + nl
            
            
            D = {}
            
            for j in range(self.numCoef[i]):
                
                if self.coef[i][j] == 0.: continue
                
                n = self.decompose_j(j, self.Nvec[i])
                
                keys = [v+str(w) for v,w in zip(namesX, n)]
                
                zm.misc.nestedDictAssign(D, keys, self.coef[i][j])
                
            
            reducedNvec = [j for j in self.Nvec[i][:-1]]
            
            n0 = [None] * (self.db.numIndVar-1)
            
            # keys = []
            # p = 0
            
            for j in range(self.calcNumCoef(reducedNvec)):
                
                n = self.decompose_j(j, reducedNvec)
                
                # print(n0, n)
                
                # input()
                
                keys = [v+str(w) for v,w in zip(namesX[:-1], n)]
                d = zm.misc.nestedDictGet(D, keys)
                keys = [k for k in d]
                
                if len(keys) == 0:
                    # n0 = n[:]
                    continue
                elif n0[0] != None:
                    # m = sum([1 for k in range(len(n0)-1) if n0[k] != n[k]])
                    
                    for k in range(len(n0)):
                        if n0[k] != n[k]:
                            break
                    m = len(n0) - k
                    
                    s += ')'*m + ' +' + nl
                
                t = tab*2
                
                flag = False
                for k in range(self.db.numIndVar-1):
                    
                    if flag or n[k] != n0[k]:
                        if n[k] != 0:
                            t += '{}{} * ('.format(namesX[k], n[k])
                        else:
                            t += '{}('.format(' '*(len(namesX[k])+4))
                        # p += 1
                        flag = True
                    else:
                        t += ' '*(len(namesX[k]) + len(str(n[k])) + 4)
                
                if keys[0][-1] != '0':
                    s += t + keys[0] + ' * {}'.format(d[keys[0]])
                else:
                    s += t + ' '*len(keys[0]) + '   {}'.format(d[keys[0]])
                
                t = ' '*len(t)
                for key in keys[1:]:
                    s += ' +' + nl
                    s += t + key + ' * {}'.format(d[key])
                
                # s += ')'
                
                # if len(keys) > 0:
                    # s += ') + ' + nl
                
                
                n0 = n[:]
            
            s += ')'*self.db.numIndVar + nl
            
        
        if not latex: return s
        
        ##############################################################
        ##############################################################
        
        S = s.split(sep='def evaluate_')[1:]
        
        # myVars = (r'\alpha', r'\beta', r'\delta_a')
        
        lS = len(S)
        nameY = [None]*lS
        nameX = [None]*lS
        eqn = [None]*lS
        front = [None]*lS
        out = [None]*lS
        
        for i in range(lS):
            s = S[i]
            for j in range(len(s)):
                if s[j] == '(':
                    J = j
                elif s[j] == ')':
                    JJ = j
                elif s[j:j+9] == 'return (\n':
                    break
            nameY[i] = s[:J]
            nameX[i] = s[J+1:j].split(sep='1, ')[:-1]
            if i < lS -1:
                eqn[i] = s[j+9:-4]
            else:
                eqn[i] = s[j+9:-2]
            
            j = -1
            for old,new in zip(nameX[i],namesX):
                j += 1
                eqn[i] = eqn[i].replace(old, new+'^')
                eqn[i] = eqn[i].replace('^1', '')
            
            front[i] = r'\quad ' + namesY[i] + r'({}'.format(namesX[0])
            for j in range(1,len(namesX)):
                front[i] += r',{}'.format(namesX[j])
            front[i] += r')=' + '\n'
            
            s = '$' + front[i] + eqn[i] + '$'
            
            s = s.replace(' '*4, r'\quad')
            
            
            s = r'\\|\quad '.join(s.split('\n'))
            
            out[i] = s
        
        return out
    
    
    
    
    












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
    
    d3 = {'maxOrder': 2, }
    
    kw = [d1, d2, d3, d1, d2, d3, d2, d1, d2, d1]
    # kw = [d3, d3, d3, d3, d3, d3, d3, d3, d3, d3]
    
    myFits = polyFit(db, kw, mpFits=1)
    
    # for c in myFits.coef:
        # print()
        # print(c)
