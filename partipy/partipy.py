import numpy as np
import flopy
import flopy.utils.binaryfile as bf


class track_particles():
	def __init__(self,cbbfile,model_ws,namfile,starting_locs,n,delt,ntimes):
		self.cbbobj = bf.CellBudgetFile(os.path.join(model_ws,modelname+'.cbc'))
		self.times = cbbobj.get_times()
		self.mf = flopy.modflow.Modflow.load(os.path.join(model_ws,namfile))
		self.delr = self.mf.dis.delr.array
		self.delc = self.mf.dis.delc.array
		self.nlay = self.mf.dis.nlay
		self.nrow = self.mf.dis.nrow
		self.ncol = self.mf.dis.ncol

	def what_cell_am_i_in(pt,nlay,nrow,ncol,delr,delc):
	    delr = delr.cumsum()
	    delc = delc.cumsum()
	    x,y = pt

	    row = np.searchsorted(delr,y)
	    col = np.searchsorted(delc,x)
	    lay = 0 # we'll get there

	    return lay, row, col

	def get_flow_faces(time,lay=0):
		idx = np.searchsorted(self.times,time)
		fff = cbb.get_data(text='FLOW FRONT FACE', totim=self.times[idx])[0]
		frf = cbb.get_data(text='FLOW RIGHT FACE', totim=slef.times[idx])[0]
		return fff, frf

	def rk4(starting_locs,n=.3,delt=20,ntimes=30):
		# Runge-Kutta implicit method
		end_pts = {}
		for i in range(len(starting_locs)):
			x0, y0 = starting_locs[i]
			print(x0,y0)
			xpts, ypts = [x0], [y0]
			for t in range(ntimes):
				xp0, yp0 = xpts[-1], ypts[-1]

				l,r,c = what_cell_am_i_in((xp0,yp0),1,nrow,ncol,delc,delr)
				fff,frf = get_flow_faces(self.times+t*delt)

				vxp0 = (frf[l][r,c]) / (delr[c] * thk * n)
				xp1 = xp0+ vxp0 * delt/2 # x0 + (vx)0 * delt

				vyp0 = (fff[l][r,c]) / (delc[r] * thk * n)
				yp1 = yp0+ vyp0 * delt/2 # y0 + (vy)0 * delt

				l,r,c = what_cell_am_i_in((xp1,yp1),1,nrow,ncol,delc,delr)
				vxp1 = (frf[l][r,c]) / (delr[c] * thk * n)
				xp2 = xp0 + vxp1 * delt/2

				vyp1 = (fff[l][r,c]) / (delc[r] * thk * n)
				yp2 = yp0 + vyp1 * delt/2

				l,r,c = what_cell_am_i_in((xp2,yp2),1,nrow,ncol,delc,delr)
				vxp2 = (frf[l][r,c]) / (delr[c] * thk * n)
				xp3 = xp0 + vxp2 * delt

				vyp2 = (fff[l][r,c]) / (delc[r] * thk * n)
				yp3 = yp0 + vyp2 * delt

			xpts.append(xp3)
			ypts.append(yp3)
			end_pts[i] = [xpts,ypts]

			return end_pts
