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

f = open(fn+'autoFit_CL.json', 'r')
clDict = json.load(f)
f.close()
f = open(fn+'autoFit_CD.json', 'r')
cdDict = json.load(f)
f.close()
f = open(fn+'autoFit_Cm.json', 'r')
cmDict = json.load(f)
f.close()

aCL, nvecCL, r2CL = pf.dict2list(clDict)
aCD, nvecCD, r2CD = pf.dict2list(cdDict)
aCm, nvecCm, r2Cm = pf.dict2list(cmDict)

f = open(fn+'auto_a5dp10.txt', 'r')
clval = float(f.readline())
cdval = float(f.readline())
cmval = float(f.readline())
f.close()

aoa, dp = 5.*np.pi/180., 10.*np.pi/180.

def test_zero_tol():
    
    aaCL, nnvecCL, rr2CL = pf.autoPolyFit(X, cl, tol=0., verbose=False)
    
    assert len(nvecCL) == len(nnvecCL)
    for v in range(len(nvecCL)):
        assert nvecCL[v] == nnvecCL[v]
    assert len(aCL) == len(aaCL)
    for j in range(pf.calcJ(nvecCL)):
        assert aCL[j] == aaCL[j]
    assert r2CL == rr2CL
    
    cclval = pf.multivariablePolynomialFunction(aCL, nvecCL, [aoa, dp])
    assert clval == cclval

def test_sigmaMultiplier():
    
    aaCD, nnvecCD, rr2CD = pf.autoPolyFit(X, cd, sigmaMultiplier=0.01, verbose=False)
    
    assert len(nvecCD) == len(nnvecCD)
    for v in range(len(nvecCD)):
        assert nvecCD[v] == nnvecCD[v]
    assert len(aCD) == len(aaCD)
    for j in range(pf.calcJ(nvecCD)):
        assert aCD[j] == aaCD[j]
    assert r2CD == rr2CD
    
    ccdval = pf.multivariablePolynomialFunction(aCD, nvecCD, [aoa, dp])
    assert cdval == ccdval

def test_MaxOrder():
    
    aaCm, nnvecCm, rr2Cm = pf.autoPolyFit(X, cm, MaxOrder=7, verbose=False)
    
    assert len(nvecCm) == len(nnvecCm)
    for v in range(len(nvecCm)):
        assert nvecCm[v] == nnvecCm[v]
    assert len(aCm) == len(aaCm)
    for j in range(pf.calcJ(nvecCm)):
        assert aCm[j] == aaCm[j]
    assert r2Cm == rr2Cm
    
    ccmval = pf.multivariablePolynomialFunction(aCm, nvecCm, [aoa, dp])
    
    assert cmval == ccmval
