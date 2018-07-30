import numpy as np
import flopy
import flopy.utils.binaryfile as bf
import os
import matplotlib.pyplot as plt

def what_cell_am_i_in(pt, nlay, delr, delc):
    delr = delr.cumsum()
    nrow = len(delc)
    delc = delc.cumsum()
    x, y = pt

    row = np.searchsorted(delr, y)
    row = nrow - np.searchsorted(delr, y) - 1  # because math
    col = np.searchsorted(delc, x)
    lay = 0  # we'll get there

    return lay, row, col

def get_flow_faces(times,time,cbbobj,lay=0):
    idx = np.searchsorted(times,time)
    if idx >= len(times):
        idx = len(times)-1
    # print(times,time,idx)
    fff = cbbobj.get_data(text='FLOW FRONT FACE', totim=times[idx])[0]
    frf = cbbobj.get_data(text='FLOW RIGHT FACE', totim=times[idx])[0]
    return fff, frf # reverse flow up and down for flow front face

def weak_sink(times,time,cbbobj,fff,frf,l,r,c):
    nlay,nrow,ncol = fff.shape
    if (r < 0) or (r >= nrow):
        return True
    if (c < 0) or (c >= ncol):
        print('batman')
        return True

    Qin = -fff[l][r,c] + frf[l][r,c] # + flf? flow from other layer, still getting there.
    idx = np.searchsorted(times,time)
    if idx >= len(times):
        idx = 0
    chb = np.asarray(cbbobj.get_data(text='Constant Head',totim=times[idx], full3D=True))[0]
    chb[chb == 0] = np.nan
    Qsnk = chb[l][r,c]
    snk = Qsnk/Qin

    if snk < 1:
        return True
    else:
        return False


class track_particles():
    def __init__(self,cbbfile,model_ws,modelname,starting_locs,n,delt,ntimes):
        self.cbbobj = bf.CellBudgetFile(os.path.join(model_ws,cbbfile))
        self.times = self.cbbobj.get_times()
        self.mf = flopy.modflow.Modflow.load(os.path.join(model_ws,modelname+'.nam'))
        self.delr = self.mf.dis.delr.array
        self.delc = self.mf.dis.delc.array
        self.nlay = self.mf.dis.nlay
        self.nrow = self.mf.dis.nrow
        self.ncol = self.mf.dis.ncol
        self.starting_locs = starting_locs
        self.n = n
        self.delt = delt
        self.ntimes = ntimes


    def rk4(self):
        # Runge-Kutta implicit method
        delr,delc = self.delr, self.delc
        thk = self.mf.dis.thickness.array
        n = self.n
        starting_locs = self.starting_locs
        delt = self.delt
        end_pts = {}
        for i in range(len(starting_locs)):
            x0, y0 = starting_locs[i]
            xpts, ypts = [x0], [y0]
            for t in range(self.ntimes):
                time = delt*t + self.times[0]
                xp0, yp0 = xpts[-1], ypts[-1]

                l,r,c = what_cell_am_i_in((xp0,yp0),1,delc,delr)
                fff,frf = get_flow_faces(self.times,time,self.cbbobj)

                vxp0 = (frf[l][r,c]) / (delr[c] * thk[l][r,c] * n)
                xp1 = xp0+ vxp0 * delt/2

                vyp0 = -(fff[l][r,c]) / (delc[r] * thk[l][r,c] * n)
                yp1 = yp0+ vyp0 * delt/2

                l,r,c = what_cell_am_i_in((xp1,yp1),1,delc,delr)
                vxp1 = (frf[l][r,c]) / (delr[c] * thk[l][r,c] * n)
                xp2 = xp0 + vxp1 * delt/2

                vyp1 = -(fff[l][r,c]) / (delc[r] * thk[l][r,c] * n)
                yp2 = yp0 + vyp1 * delt/2

                l,r,c = what_cell_am_i_in((xp2,yp2),1,delc,delr)
                if (r < 0) or (r >= fff.shape[1]):
                    break
                if (c < 0) or (c >= fff.shape[2]):
                    break
                vxp2 = (frf[l][r,c]) / (delr[c] * thk[l][r,c] * n)
                xp3 = xp0 + vxp2 * delt


                vyp2 = -(fff[l][r,c]) / (delc[r] * thk[l][r,c] * n)
                yp3 = yp0 + vyp2 * delt

                l,r,c, = what_cell_am_i_in((xp3,yp3),1,delc,delr)
                if weak_sink(self.times,time,self.cbbobj,fff,frf,l,r,c):
                    break
                else:
                    xpts.append(xp3)
                    ypts.append(yp3)
            end_pts[i] = [xpts,ypts]

        return end_pts
