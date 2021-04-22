import matplotlib.pyplot as plt 
import numpy as np
import h5py

def plot_spectrum(ds,cut):
    qs=ds['q'][cut]
    nbins=int(np.max(qs)/1000)
    hist,bins=np.histogram(qs,bins=nbins)
    bins=0.5*(bins[1:] + bins[:-1])
    fig=plt.figure()
    plt.step(bins,hist,where='mid')
    plt.ylabel('#')
    plt.xlabel('q')
    plt.yscale('log')
    plt.grid()
    plt.show()
    return 0

def plot_correlation(ds,x,y,cut):
    fig=plt.figure()
    plt.plot(ds[x][cut],ds[y][cut],'.')
    plt.xlabel(x)
    plt.ylabel(y)
    plt.grid()
    plt.show()
    return 0


# ========================================
# http://docs.h5py.org/en/stable/quick.html
# https://www.pythonforthelab.com/blog/how-to-use-hdf5-files-in-python/
# https://www.christopherlovell.co.uk/blog/2016/04/27/h5py-intro.html

# input file:
fevents="/Users/kiwi/Documents/PROYECTOS/tesis_nicolas_ib/events.hdf5"
f = h5py.File(fevents,"r")
print("Groups in the HDF5: ",f.keys())
print("Datasets in group 'hits': ",f['hits'].keys())

# Print datasets in the HDF5:
def printname(name):
    print(name)
f.visit(printname) 

# Load a dataset
ds=f['hits/parameters'] # ds becomes a numpy array. 
print("dataset shape: ",ds.shape) # Print number of events
print("dataset type: ",ds.dtype) 
#print(ds[...]) # Print all dataset content
#print(ds['q'][...]) # Print only the 'q' column of the dataset

# Spectrum plot
cut=(ds['q']>5.) & (ds['ohdu']==0) & (ds['flag']==0) 
plot_spectrum(ds,cut)

# Correlation plot
cut=(ds['q']<1.8e6) & (ds['ohdu']==0) & (ds['flag']==0) 
plot_correlation(ds,'q','npix',cut)

## Create a HDF5 with a subset:
#cut=(ds['q']<1.8e6) & (ds['ohdu']==0) & (ds['flag']==0) 
#fnew = h5py.File("events_clean.hdf5", "w")
#f.copy('imgs_metadata',fnew)
#hits_group=fnew.create_group('hits')
#for ds_name in f['hits'].keys(): 
#    print(ds_name)
#    hits_group.create_dataset(i,f["hits/"+ds_name][cut],dtype=f[ds_name].dtype)
#fnew.close()
##
exit()

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


