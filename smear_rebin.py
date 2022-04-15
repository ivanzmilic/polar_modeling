import sys
import matplotlib.pyplot as plt 
from astropy.io import fits 
import numpy as np 

from scipy.ndimage import gaussian_filter

file_in = sys.argv[1]
file_out = sys.argv[2]

D = float(sys.argv[3])
binning = int(sys.argv[4])

data_in = fits.open(file_in)

test = data_in[0].data[:,:,0,0]
print (test.shape)

NL,NX= test.shape

test = data_in[0].data[0,0,:,:]

NY,NS= test.shape

print (NX,NY, NS, NL)

NX_new = NX // binning
NY_new = NY // int(binning / np.cos(np.radians(65)))

print (NX_new,NY_new)

stokes_out = np.zeros([NX_new,NY_new,NS,NL])

sigma = 1.22 * 6300E-10 / D * 206265 * 720. / 14. / 2.35
print ("PSF in pixels is: ", sigma)

for s in range (0,NS):
	for l in range(0,NL):

		frame = np.copy(data_in[0].data[l,:,:,s])

		frame = gaussian_filter(frame,(sigma,sigma/np.cos(np.radians(65))))

		# Bin

		binx = binning
		biny = int(binning / np.cos(np.radians(65)))

		for i in range(0,NX_new):
			for j in range(0,NY_new):
				stokes_out[i,j,s,l] = np.mean(frame[i*binx:(i+1)*binx,j*biny:(j+1)*biny])

myhdu = fits.PrimaryHDU(stokes_out)
myhdu.writeto(file_out,overwrite=True)


