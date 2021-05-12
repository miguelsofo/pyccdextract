import matplotlib.pyplot as plt 
import matplotlib as matplt
import numpy as np
import h5py

def plot_spectrum(ds,cut):
    qs=ds['q'][cut]
    nbins=int(qs.size/100)
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

def plot_events(f,events):
    # Plot the events
    xpix=f['hits/xpix']
    ypix=f['hits/ypix']
    pix=f['hits/pix']
    xpix=xpix[events]
    ypix=ypix[events]
    pix=pix[events]
    di=f['imgs_metadata']
    img=np.zeros((di['nrow'][0],di['ncol'][0]))
    img[np.concatenate(ypix),np.concatenate(xpix)]=np.concatenate(pix)
    
    fig=plt.figure()
    cmap=matplt.cm.jet
    cmap.set_under('white')
    plt.imshow(img,vmin=0.1,vmax=np.max(img),cmap=cmap)
    cbar=plt.colorbar()
    cbar.set_label('value',fontsize=14)
    plt.grid()
    plt.xlabel('x (pixel)',fontsize=14)
    plt.ylabel('y (pixel)',fontsize=14)
    plt.show()
    return 0

# ========================================
# http://docs.h5py.org/en/stable/quick.html
# https://www.pythonforthelab.com/blog/how-to-use-hdf5-files-in-python/
# https://www.christopherlovell.co.uk/blog/2016/04/27/h5py-intro.html

# input file:
fevents="/Volumes/ExtremeSSD/datos-ACDS/test/event.hdf5"

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

ohdu=0
# Spectrum plot
cut=(ds['q']<5e7) & (ds['ohdu']==ohdu) & (ds['flag']==0) 
plot_spectrum(ds,cut)

# Correlation plot
cut=(ds['q']<1.8e6) & (ds['ohdu']==ohdu) & (ds['flag']==0) 
plot_correlation(ds,'q','npix',cut)

plot_events(f,cut)

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


