## Interpolate SPECFEM seismograms to new sampling interval
#
#  Daniel Koehn
#  Kiel, 22.01.2018

import numpy as np
from scipy import interpolate
 
# read seismograms 
basename_in = "seis_MXZ/"
basename_out = "seis_MXZ_dt_1s/"
state = "EU"

minlat = 35.5
maxlat = 58.5

minlon = -14
maxlon = 30

# define longitude and latitude positions of receivers
drec = 0.45   # receiver distance [degree]
lat_rec = np.arange(minlat, maxlat, drec)
lon_rec = np.arange(minlon, maxlon, drec)

for i in range(0,lat_rec.size,1):
	for j in range(0,lon_rec.size,1):
	
		station = "H" + "%0.*f" %(0,np.fix(i)) + "R" + "%0.*f" %(0,np.fix(j))
		
		# define file names for z-component
		name_seis_Z = basename_in + state + "." + station + ".MXZ.sem.ascii"
    
		# read z-component data
		data_Z = np.genfromtxt(name_seis_Z, dtype=None)

		# define time vector [s]
		t = data_Z[:,0]

		# Interpolate seismograms to 1 s
		dt = 1.0   # new sampling interval [s]
		tnew = np.arange(np.min(t),np.max(t),dt)
		ns = np.size(tnew)
		tmp = data_Z[:,1]
		f = interpolate.interp1d(t, tmp)
		seis_new = f(tnew)    

		# output interpolated seismograms
		name_seis_out_Z = basename_out + state + "." + station + ".MXZ.sem.ascii"

		# assemble vectors into an array
		tmp = np.zeros(ns, dtype=[('var1', float), ('var2', float)])
		tmp['var1'] = tnew
		tmp['var2'] = seis_new

		# write seismogram to file
		np.savetxt(name_seis_out_Z, tmp, fmt='%5.3f %5.3e')
