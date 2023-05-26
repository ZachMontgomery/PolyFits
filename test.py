import numpy as np
import polyFits as pf
import ZachsModules as zm
zm.zp.updateRCParams()
plt = zm.plt

N = 3
x1 = np.linspace(-10,10,N)

y = x1 ** 2
db = pf.database(x1, y)

N = 100
tol = 1e-6
x = np.linspace(min(x1)+tol,max(x1)-tol,N)
y = []
prog = zm.io.oneLineProgress(N)
for z in x:
    y.append(db.interpolate([z]))
    prog.display()

plt.plot(db.x, db.y, 'o')
plt.plot(x, y, '-k')

plt.show()

# db = pf.database([[0], [1]], [[2],[3]])

# print(db.interpolate([0.5]))
