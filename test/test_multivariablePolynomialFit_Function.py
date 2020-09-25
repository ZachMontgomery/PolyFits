import numpy as np
import polyFits as pf
import json

fn = './test/'
f = open(fn+'database.txt', 'r')
database = f.readlines()
f.close()

aoa, dp, cl, cd, cm = [], [], [], [], []
for line in database[1:]:
    aoa.append( float( line[  8: 25] ) )
    dp.append(  float( line[ 34: 51] ) )
    cl.append(  float( line[ 60: 77] ) )
    cd.append(  float( line[ 87:103] ) )
    cm.append(  float( line[112:   ] ) )
X = np.array([aoa, dp]).T

f = open(fn+'fit_CL.json', 'r')
clDict = json.load(f)
f.close()
f = open(fn+'fit_CD.json', 'r')
cdDict = json.load(f)
f.close()
f = open(fn+'fit_Cm.json', 'r')
cmDict = json.load(f)
f.close()

aCL, nvecCL, r2CL = pf.dict2list(clDict)
aCD, nvecCD, r2CD = pf.dict2list(cdDict)
aCm, nvecCm, r2Cm = pf.dict2list(cmDict)

f = open(fn+'a5dp10.txt', 'r')
clval = float(f.readline())
cdval = float(f.readline())
cmval = float(f.readline())
f.close()

aoa, dp = 5.*np.pi/180., 10.*np.pi/180.

def test_simpleConstriants():
    
    aaCL, rr2CL = pf.multivariablePolynomialFit(nvecCL, X, cl, sym_same=[(0,1)], verbose=False)
    
    assert len(aCL) == len(aaCL)
    for j in range(pf.calcJ(nvecCL)):
        assert aCL[j] == aaCL[j]
    assert r2CL == rr2CL
    
    cclval = pf.multivariablePolynomialFunction(aCL, nvecCL, [aoa, dp])
    assert clval == cclval

def test_percent():
    
    aaCD, rr2CD = pf.multivariablePolynomialFit(nvecCD, X, cd, sym_diff=[(0,1)], percent=True, verbose=False)
    
    assert len(aCD) == len(aaCD)
    for j in range(pf.calcJ(nvecCD)):
        assert aCD[j] == aaCD[j]
    assert r2CD == rr2CD
    
    ccdval = pf.multivariablePolynomialFunction(aCD, nvecCD, [aoa, dp])
    assert cdval == ccdval

def test_weighting():
    
    def w(x, y, p):
        if abs(y[p]) < 0.0001:
            return 1.
        return 0.0001 / abs(y[p])
    
    aaCm, rr2Cm = pf.multivariablePolynomialFit(nvecCm, X, cm, sym_same=[(0,1)], weighting=w, verbose=False)
    
    assert len(aCm) == len(aaCm)
    for j in range(pf.calcJ(nvecCm)):
        assert aCm[j] == aaCm[j]
    assert r2Cm == rr2Cm
    
    ccmval = pf.multivariablePolynomialFunction(aCm, nvecCm, [aoa, dp])
    assert cmval == ccmval

