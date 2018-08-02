import flopy
import pandas as pd
import numpy as np
import os


import numpy as np
import flopy
import os
import platform
# Assign name and create modflow model object
modelname = 'tutorial1'
exe = os.path.join('gw_codes','mf2005.exe')
if platform.system() == 'Darwin':
    exe = 'mf2005' # assuming you have mfnwt in your path
model_ws = os.path.join('presentation')
mf = flopy.modflow.Modflow(modelname, exe_name=exe,model_ws=model_ws)

# Model domain and grid definition
Lx = 1000.
Ly = 1000.
ztop = 0.
zbot = -50.
nlay = 1
nrow = 11
ncol = 11
delr = Lx/ncol
delc = Ly/nrow
delv = (ztop - zbot) / nlay
botm = np.linspace(ztop, zbot, nlay + 1)
nper=3
# Create the discretization object
dis = flopy.modflow.ModflowDis(mf, nlay, nrow, ncol, delr=delr, delc=delc,
                               top=ztop, botm=botm[1:],nper=nper,perlen=365)

# Variables for the BAS package
ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
ibound[:, 0, :] = -1
ibound[:, -1, :] = -1
strt = np.ones((nlay, nrow, ncol), dtype=np.float32)
strt[:, 0, :] = 10.
strt[:, -1, :] = 9.5
bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt)

# Add LPF package to the MODFLOW model
lpf = flopy.modflow.ModflowLpf(mf, hk=50., vka=10., ipakcb=53)

# Add OC package to the MODFLOW model
spd = {}
for sp in range(nper):
    spd[(sp, 0)] = ['print head', 'print budget', 'save head', 'save budget']
oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd, compact=True)

# Add PCG package to the MODFLOW model
pcg = flopy.modflow.ModflowPcg(mf)

wel = flopy.modflow.ModflowWel(mf,stress_period_data={1:[0,5,5,1000],2:[0,5,5,-10000]})

# Write the MODFLOW model input files
mf.write_input()

# Run the MODFLOW model
success, buff = mf.run_model()

# Post process the results
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf

plt.subplot(1, 1, 1, aspect='equal')
hds = bf.HeadFile(os.path.join(model_ws,modelname + '.hds'))
times = hds.get_times()
head = hds.get_data(totim=times[-1])
levels = np.arange(1, 10, 1)
extent = (delr / 2., Lx - delr / 2., Ly - delc / 2., delc / 2.)
plt.contour(head[0, :, :], levels=levels, extent=extent)
plt.savefig(os.path.join(model_ws,'tutorial1a.png'))

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1, 1, 1, aspect='equal')

hds = bf.HeadFile(os.path.join(model_ws,modelname+'.hds'))
head = hds.get_data(totim=times[0])
levels = np.linspace(0, 10, 11)

cbb = bf.CellBudgetFile(os.path.join(model_ws,modelname+'.cbc'))
kstpkper_list = cbb.get_kstpkper()
frf = cbb.get_data(text='FLOW RIGHT FACE', totim=times[0])[0]
fff = cbb.get_data(text='FLOW FRONT FACE', totim=times[0])[0]

modelmap = flopy.plot.ModelMap(model=mf, layer=0)
qm = modelmap.plot_ibound()
lc = modelmap.plot_grid()
# cs = modelmap.contour_array(head, levels=levels)
quiver = modelmap.plot_discharge(frf, fff, head=head)
plt.title('Flow Before Well',fontsize=25)
fig.tight_layout()
fig.savefig(os.path.join(model_ws,'Steady_state.png'))

##########################################
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1, 1, 1, aspect='equal')

hds = bf.HeadFile(os.path.join(model_ws,modelname+'.hds'))
head = hds.get_data(totim=times[1])
levels = np.linspace(0, 10, 11)

cbb = bf.CellBudgetFile(os.path.join(model_ws,modelname+'.cbc'))
kstpkper_list = cbb.get_kstpkper()
frf = cbb.get_data(text='FLOW RIGHT FACE', totim=times[1])[0]
fff = cbb.get_data(text='FLOW FRONT FACE', totim=times[1])[0]

modelmap = flopy.plot.ModelMap(model=mf, layer=0)
qm = modelmap.plot_ibound()
lc = modelmap.plot_grid()
wc = modelmap.plot_bc('wel',color='r',kper=2)

quiver = modelmap.plot_discharge(frf, fff, head=head)
plt.title('Flow with Injection Well',fontsize=25)
fig.tight_layout()
fig.savefig(os.path.join(model_ws,'injection.png'))

##########################################
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(1, 1, 1, aspect='equal')

hds = bf.HeadFile(os.path.join(model_ws,modelname+'.hds'))
head = hds.get_data(totim=times[2])
levels = np.linspace(0, 10, 11)

cbb = bf.CellBudgetFile(os.path.join(model_ws,modelname+'.cbc'))
kstpkper_list = cbb.get_kstpkper()
frf = cbb.get_data(text='FLOW RIGHT FACE', totim=times[2])[0]
fff = cbb.get_data(text='FLOW FRONT FACE', totim=times[2])[0]

modelmap = flopy.plot.ModelMap(model=mf, layer=0)
qm = modelmap.plot_ibound()
lc = modelmap.plot_grid()
wc = modelmap.plot_bc('wel',color='r',kper=2)
# cs = modelmap.contour_array(head, levels=levels)
quiver = modelmap.plot_discharge(frf, fff, head=head)
plt.title('Flow with Recovery Well',fontsize=25)
fig.tight_layout()
fig.savefig(os.path.join(model_ws,'extraxtion.png'))

plt.show()