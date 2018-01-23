## DFT of SPECFEM seismograms for one frequency
#
#  Daniel Koehn
#  Kiel, 22.01.2018

import numpy as np
 
# read seismograms 
basename_in = "OUTPUT_FILES/seis_MXZ_dt_1s/"
state = "EU"

minlat = 35.5
maxlat = 58.5

minlon = -14
maxlon = 30

# define longitude and latitude positions of receivers
drec = 0.45   # receiver distance [degree]
lat_rec = np.arange(minlat, maxlat, drec)
lon_rec = np.arange(minlon, maxlon, drec)

# array of DFT for each station
DFT_real = np.ones((lat_rec.size,lon_rec.size))

# loop over all stations
for i in range(0,lat_rec.size,1):
	for j in range(0,lon_rec.size,1):
	
		station = "H" + "%0.*f" %(0,np.fix(i)) + "R" + "%0.*f" %(0,np.fix(j))
		
		# define file names for z-component
		name_seis_Z = basename_in + state + "." + station + ".MXZ.sem.ascii"
    
		# read z-component data
		data_Z = np.genfromtxt(name_seis_Z, dtype=None)

		# define time vector t [s] and number of samples ns
		t = data_Z[:,0]
		ns = np.size(t)
		
		# time sampling DT [s]
		DT = t[1] - t[0]
		
		# FD at frequency f [Hz] to extract from seismogram
		f = 1.0 / 20.0

		DFT_real[i][j] = 0.0

		# Calculate DFT of seismogram
		for k in range(0,ns,1):			
			trig = np.cos(2.0 * np.pi * t[k] * f) * DT
			DFT_real[i][j] += trig * data_Z[k,1]
			
# output of DFT data
file_out = "DFT_real.bin"
f = open (file_out, mode='wb')
data_type = np.dtype ('float32').newbyteorder ('<')
DFT_real = np.array(DFT_real, dtype=data_type)
DFT_real = np.rot90(DFT_real,1)
DFT_real = np.flipud(DFT_real)
DFT_real.tofile(f)
f.close()
