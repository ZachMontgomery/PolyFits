if __name__ == '__main__':
    import numpy as np
    import polyFits as pf
    import ZachsModules as zm
    from mpl_toolkits.mplot3d import axes3d
    
    def evalMultiPoly(A, x, Nvec):
        a = A[:]
        for v in range(len(Nvec)-1,0,-1):
            k = Nvec[v] + 1
            J = pf.polyFit.calcNumCoef(Nvec[:v+1])
            for i in range(int(J/k)):
                s = i * k
                e = s + k
                a[i] = pf.polyFit.evalPoly1D(a[s:e], x[v])
        return pf.polyFit.evalPoly1D(a[:int(J/k)], x[0])
    
    Nvec = [2, 3, 4]
    V = len(Nvec)
    J = pf.polyFit.calcNumCoef(Nvec)
    a = [1.0] * J
    n = 11
    erMag = 1.0
    numPoints = n ** V
    x = np.zeros((numPoints, V))
    y = np.zeros((numPoints, 4))
    error = (np.random.rand(numPoints) - 0.5) * erMag
    i = -1
    prog = zm.io.oneLineProgress(numPoints, msg='Creating dummy data')
    for xs in pf.nestedFor(*[np.linspace(-1,1,n)]*V):
        i += 1
        y[i,:2] = evalMultiPoly(a, xs, Nvec)
        y[i,2:] = evalMultiPoly(a, xs, Nvec) + error[i]
        x[i,:] = xs
        prog.display()
    
    db = pf.database(x, y, namesY=('manual fit, no error', 'auto fit, no error', 'manual fit, error', 'auto fit, error'))
    
    zm.zp.updateRCParams(ax3Dsetup=True)
    plt = zm.plt
    plt.ion()
    
    fig, ax = plt.subplots(2,2, subplot_kw={'projection':'3d'}, figsize=(6.4*2, 4.8*2))
    ax = np.reshape(ax, -1)
    for a in ax: zm.zp.axis3DgridZachsPreferences(a)
    zm.zp.link3Daxes(fig, ax)
    
    db.viewData(fig, ax, (0,1,2,3), )
    
    d1 = {'Nvec': Nvec}
    d2 = {}
    kw = [d1, d2, d1, d2]
    
    fit = pf.polyFit(db, kw, )
    
    print(fit.writeHardCodedEqs())
    
    fig, ax = plt.subplots(2,2, subplot_kw={'projection':'3d'}, figsize=(6.4*2, 4.8*2))
    ax = np.reshape(ax, -1)
    for a in ax: zm.zp.axis3DgridZachsPreferences(a)
    zm.zp.link3Daxes(fig, ax)
    
    db.viewData(fig, ax, (0,1,2,3), f=fit.f)
    
