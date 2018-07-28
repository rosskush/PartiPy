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
delr = Lx/ncol * np.ones(ncol)
delc = Ly/nrow * np.ones(nrow)
delv = (ztop - zbot) / nlay
botm = np.linspace(ztop, zbot, nlay + 1)
nper=3
thk = ztop - zbot
def what_cell_am_i_in(pt,nlay,nrow,ncol,delr,delc):
    delr = delr.cumsum()
    delc = delc.cumsum()
    x,y = pt

    row = np.searchsorted(delr,y)
    col = np.searchsorted(delc,x)
    lay = 0

    return lay, row, col

def track_particle(starting_locs,n=.3,delt=30,ntimes=50):
    end_pts = {}
    for i in range(len(starting_locs)):
        x0, y0 = starting_locs[i]
        xpts, ypts = [x0], [y0]
        for t in range(ntimes):
            xp0, yp0 = xpts[-1], ypts[-1]

            l,r,c = what_cell_am_i_in((xp0,yp0),1,nrow,ncol,delc,delr)
            print(xp0,yp0)
            print(l,r,c)
            vxp0 = (frf[l][r,c]) / (delr[c] * thk * n)
            xp1 = xp0 + vxp0 * delt/2 # x0 + (vx)0 * delt

            vyp0 = (fff[l][r,c]) / (delc[r] * thk * n)
            yp1 = yp0 - vyp0 * delt/2 # y0 + (vy)0 * delt


            l,r,c = what_cell_am_i_in((xp1,yp1),1,nrow,ncol,delc,delr)
            print(xp1,yp1)
            print(l,r,c)
            vxp1 = (frf[l][r,c]) / (delr[c] * thk * n)
            xp2 = xp0 + vxp1 * delt/2

            vyp1 = (fff[l][r,c]) / (delc[r] * thk * n)
            yp2 = yp0 -vyp1 * delt/2


            l,r,c = what_cell_am_i_in((xp2,yp2),1,nrow,ncol,delc,delr)
            print(xp2,yp2)
            print(l,r,c)
            vxp2 = (frf[l][r,c]) / (delr[c] * thk * n)
            xp3 = xp0 + vxp2 * delt

            vyp2 = (fff[l][r,c]) / (delc[r] * thk * n)
            yp3 = yp0 - vyp2 * delt

            l,r,c = what_cell_am_i_in((xp3,yp3),1,nrow,ncol,delc,delr)
            print(xp3,yp3)
            print(l,r,c)
            # exit()

            xpts.append(xp3)
            ypts.append(yp3)
        end_pts[i] = [xpts,ypts]

    return end_pts

# for x
	# for t in [0,]
starting_locs = [(199,250),(150,850),(400,510)]

particles = track_particle(starting_locs)
px,py = particles[0]
px1,py1 = particles[1]
px2,py2 = particles[2]

# print(px,py)




# fig, ax = plt.subplots()
# plt.imshow(frf[0])
# plt.colorbar()
# ax.scatter(starting_pt[0],starting_pt[1])
# plt.title('FLOW RIGHT FACE')

# fig, ax = plt.subplots()
# plt.imshow(fff[0])
# plt.colorbar()
# ax.scatter(starting_pt[0],starting_pt[1])
# plt.title('FLOW FRONT FACE')

mf = flopy.modflow.Modflow.load(os.path.join(model_ws,'tutorial1.nam'))
fig, ax= plt.subplots()
modelmap = flopy.plot.ModelMap(model=mf, layer=0)
qm = modelmap.plot_ibound()
lc = modelmap.plot_grid()
hds = bf.HeadFile(os.path.join(model_ws, modelname + '.hds'))
head = hds.get_data(totim=30)
# levels = np.linspace(0, 10, 11)
# cs = modelmap.contour_array(head, levels=levels)
ax.scatter(px, py)
ax.scatter(px1,py1)
ax.scatter(px2,py2)
quiver = modelmap.plot_discharge(frf, fff, head=head)

plt.show()
