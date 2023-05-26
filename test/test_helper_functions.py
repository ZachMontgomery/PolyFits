import polyFits as pf

nvec = [6,7,8,9,10,11,12]
V = len(nvec)

def test_calcJ():
    assert pf.calcJ(nvec) == 8648640

def test_decompose_j():
    temp = pf.decompose_j(8648639, nvec)
    for i in range(V):
        assert temp[i] == nvec[i]
    temp = pf.decompose_j(0, nvec)
    for i in range(V):
        assert temp[i] == 0
    j = 4364472
    temp = pf.decompose_j(j, nvec)
    x = [3,4,2,3,4,4,8]
    for i in range(V):
        assert temp[i] == x[i]

def test_compose_j():
    temp = pf.compose_j(nvec, nvec)
    assert temp == 8648639
    temp = pf.compose_j([0]*V, nvec)
    assert temp == 0
    temp = pf.compose_j([3,4,2,3,4,4,8], nvec)
    assert temp == 4364472

def test_kDecompose():
    temp = pf.kDecompose(1254678, V)
    x = [5, 7, 1, 2, 0, 3, 6]
    for i in range(V):
        assert temp[i] == x[i]

def test_kCompose():
    temp = pf.kCompose([5, 7, 1, 2, 0, 3, 6])
    assert temp == 1254678

