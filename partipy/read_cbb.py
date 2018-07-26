import flopy
import flopy.utils.binaryfile as bf
import os
import matplotlib.pyplot as plt

import numpy as np


model_ws = os.path.join('..','examples')
modelname = 'tutorial1'

cbb = bf.CellBudgetFile(os.path.join(model_ws,modelname+'.cbc'))

times = cbb.get_times()
print(times)
kstpkper_list = cbb.get_kstpkper()
frf = cbb.get_data(text='FLOW RIGHT FACE', totim=times[-1])[0] # Flow Right Face (flow from the column to the right--higher column number),
fff = cbb.get_data(text='FLOW FRONT FACE', totim=times[-1])[0] # Flow Front Face (flow from the row below--higher row number)
# fff = cbb.get_data(text='FLOW LOWER FACE', totim=times[-1])[0] # Flow Lower Face (flow from layer below--higher layer number)
starting_pt = (1,0)

# Model domain and grid definition
Lx = 1000.
Ly = 1000.
ztop = 0.
zbot = -50.
nlay = 1
nrow = 10
ncol = 10
delr = Lx/ncol
delc = Ly/nrow
delv = (ztop - zbot) / nlay
botm = np.linspace(ztop, zbot, nlay + 1)
nper=3
thk = ztop - zbot
print(delr)
# def what_cell_am_i_in(pt,delr,delc):
	# return lay, row, col

starting_locs = [(1,0)]
def track_particle(starting_locs,n=.3,delt=.5):
    for i in range(len(starting_locs)):
        x0, y0 = starting_locs[i]
        print(x0,y0)
        # for t in times:
        vxp1 = (frf[0][0,1]) / (delr * thk * n)
        xp1 = x0+ vxp1 * delt/2 # x0 + (vx)0 * delt

        vxp2 = (frf[0][0,1]) / (delr * thk * n)
        xp2 = xp1 + vxp2 * delt/2

        vxp3 = (frf[0][0,1]) / (delr * thk * n)
        xp3 = xp2 + vxp3 * delt/2

        vxp4 = (frf[0][0,1]) / (delr * thk * n)
        xp4 = xp3 + vxp4 * delt/2


        vyp1 = (fff[0][0,1]) / (delc * thk * n)
        yp1 = y0+ vyp1 * delt/2 # y0 + (vy)0 * delt

        vyp2 = (fff[0][0,1]) / (delc * thk * n)
        yp2 = yp1 + vyp2 * delt/2

        vyp3 = (fff[0][0,1]) / (delc * thk * n)
        yp3 = yp2 + vyp3 * delt/2

        vyp4 = (fff[0][0,1]) / (delc * thk * n)
        yp4 = yp3 + vyp4 * delt/2

        print(x0, xp1, xp2, xp3, xp4)

        print(y0, yp1, yp2, yp3, yp4)

        return [x0, xp1, xp2, xp3, xp4], [y0, yp1, yp2, yp3, yp4]

# for x
	# for t in [0,]

px, py = track_particle(starting_locs)








fig, ax = plt.subplots()
plt.imshow(frf[0])
plt.colorbar()
ax.scatter(starting_pt[0],starting_pt[1])
plt.title('FLOW RIGHT FACE')

fig, ax = plt.subplots()
plt.imshow(fff[0])
plt.colorbar()
ax.scatter(starting_pt[0],starting_pt[1])
ax.scatter(px, py)
plt.title('FLOW FRONT FACE')


plt.show()
