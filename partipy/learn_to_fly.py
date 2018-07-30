import flopy
import flopy.utils.binaryfile as bf
import os
import matplotlib.pyplot as plt
import numpy as np
import partipy

model_ws = os.path.join('..','examples')
modelname = 'tutorial1'

starting_locs = [(199,250),(150,850),(400,510)]

# particles = partipy.track_particles.rk4(starting_locs=starting_locs,n=3,delt=20,ntimes=30)
tp = partipy.track_particles(modelname+'.cbc',model_ws,modelname,starting_locs,.3,20,30)

particles = tp.rk4()

px, py = particles[0]
px1,py1 = particles[1]
px2,py2 = particles[2]

mf = flopy.modflow.Modflow.load(os.path.join(model_ws,'tutorial1.nam'))
cbb = bf.CellBudgetFile(os.path.join(model_ws,modelname+'.cbc'))
times = cbb.get_times()

frf = cbb.get_data(text='FLOW RIGHT FACE', totim=times[-1])[0] # Flow Right Face (flow from the column to the right--higher column number),
fff = cbb.get_data(text='FLOW FRONT FACE', totim=times[-1])[0] # Flow Front Face (flow from the row below--higher row number)

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
plt.title('Runge-Kutta Method')

fig.savefig('partipy_example_rk.png')

plt.show()
