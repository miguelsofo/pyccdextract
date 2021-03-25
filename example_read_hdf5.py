import numpy as np
import h5py

# ========================================
# With h5py:
# (check: http://docs.h5py.org/en/stable/quick.html)
f = h5py.File("events.hdf5", "r")
# Print HDF5 organization:
def printname(name):
    print(name)
#f.visit(printname) # creo que ejecuta printname sobre todos los parametros
# Load a dataset
ds=f['hits/parameters']
#print("data set shape: ",ds.shape)
# ... From now, ds can be used as a numpy array ... 
#print(ds[...]) # Print all dataset content
#print(ds['q'][...]) # Print only the 'q' column of the dataset
#sizes=ndimage.sum(image,label_im,index=range(1,nb_labels+1))
#hist,bins=np.histogram(sizes)
#bins = 0.5*(bins[1:] + bins[:-1]) # center of the bins.

# Boolean indexing
cut=ds['q']>5. 
qs=ds['q'][cut]

# Arrays with event information
dsxpix=f['hits/xpix']
#print(dsxpix.shape)
nevent=41
q=ds['q'][nevent] # del evento 40
npix=ds['npix'][nevent]
xpix=dsxpix[nevent]
print(qs)
#print(npix)
#print(xpix)

# Plot event
#print(ds['q'][ds['q']>4.0])
#hits = f['hits']
#print(hits.shape)
#print(hits.dtype)
#print(hits[0])
f.close()


