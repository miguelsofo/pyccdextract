from scipy import ndimage
from astropy.io import fits
from datetime import datetime
import matplotlib.pyplot as plt 
import numpy as np
import math
import sys
import h5py
import os.path

def getsatvalues(imgsmtdt):
    satvalues=[]
    nhdu=len(fits.open(imgsmtdt[0][0]))
    hdumax=np.zeros(nhdu)
    print("Getting saturation values...")
    for img in imgsmtdt: #Loop over each image
        hdul = fits.open(img[0]) # read the image
        for hdu in range(0,nhdu): # loop over image hdus
            tmpmax=np.max(hdul[hdu].data)
            if tmpmax>hdumax[hdu]:
                hdumax[hdu]=tmpmax
    hdumax=hdumax*0.95
    print(hdumax)
    print("...done\n")
    return hdumax

def runExtract(config):
    # Read configuration:
    n=config["noise"] # Noise of each CCD amplifier (extension). 
    try: 
        g=config["gain"] # (ADU/e) Gain of each CCD amplifier (extension).
    except:
        print("GAIN values not provied.")
        print("Extraction of regular-CCD images.")
        g=np.ones(len(n))*np.nan
    imgsmtdt=get_images(config["in_folder"]) # Input images.
    outfile=config["out_file"] # Output file.
    # Start extraction process: 
    initOutFile(outfile) # Check and init the output HDF5 file.
    # Loop over the input images: 
    satvalues=getsatvalues(imgsmtdt)
    print("Saturation values in electrons: ")
    print(satvalues/g)
    for img in imgsmtdt:
        imgname=img[0]
        print("\n===============================================================")
        print(imgname)
        hdul = fits.open(imgname)
        if len(hdul)>len(g): # Check if the number of HDUs agrees with the number of gains.
            print("This image have more HDUs than gains in the config file.")
            print("Using only first {} ones out of {}.".format(len(g),len(hdul)));
        imgid=img[8] # image ID from image metadata.
        for hdu in range(0,len(g)):
            T1,T2=getThr(g,n,hdu)
            events=extract_hdu(imgid,hdul,hdu,T1,T2,g[hdu],satvalues[hdu])
            fillHDF5(outfile,events)
    fillHDF5_imgsmtdt(outfile,imgsmtdt) # Fill HDF5 with images metadata.
    return 0

def initHDF5(filename):
    with h5py.File(filename, "w") as f:
        # Definition of the data-types of the data-set: 
        # (check: http://docs.h5py.org/en/stable/special.html)
        # You must have h5py > v2.3, or run with python3.6.
        dt_q=np.dtype('float32')
        dt_qe=np.dtype('float32')
        dt_qq=np.dtype('float32')
        dt_npix=np.dtype('int32')
        dt_flag=np.dtype('int8')
        dt_ohdu=np.dtype('int8')
        dt_imgid=np.dtype('int32')
        dt_xmin=np.dtype('int32')
        dt_xmax=np.dtype('int32')
        dt_ymin=np.dtype('int32')
        dt_ymax=np.dtype('int32')
        dt_xbary=np.dtype('float32')
        dt_ybary=np.dtype('float32')
        dt_xvar=np.dtype('float32')
        dt_yvar=np.dtype('float32')
        dt_param=np.dtype([('q',dt_q),('qe',dt_qe),('qq',dt_qq),('npix',dt_npix),('flag',dt_flag),('ohdu',dt_ohdu),('imgid',dt_imgid),('xmin',dt_xmin),('xmax',dt_xmax),('ymin',dt_ymin),('ymax',dt_ymax),('xbary',dt_xbary),('ybary',dt_ybary),('xvar',dt_xvar),('yvar',dt_yvar)])
        # ... variable length data-types:
        dt_xpix=h5py.vlen_dtype(np.dtype('int32'))
        dt_ypix=h5py.vlen_dtype(np.dtype('int32'))
        dt_pix=h5py.vlen_dtype(np.dtype('float32'))
        dt_epix=h5py.vlen_dtype(np.dtype('float32'))
        # Create the *hits* group: it will have 4 datasets: the parameters and the xpix, ypix and pix arrays.
        hits = f.create_group('hits')
        ds_param = hits.create_dataset('parameters',(0,),dtype=dt_param,maxshape=(None,),chunks=(500,),compression="gzip",compression_opts=9,shuffle=True)
        ds_xpix = hits.create_dataset('xpix',(0,),dtype=dt_xpix,maxshape=(None,),chunks=(500,),compression="gzip",compression_opts=9,shuffle=True)
        ds_ypix = hits.create_dataset('ypix',(0,),dtype=dt_ypix,maxshape=(None,),chunks=(500,),compression="gzip",compression_opts=9,shuffle=True)
        ds_pix = hits.create_dataset('pix',(0,),dtype=dt_pix,maxshape=(None,),chunks=(500,),compression="gzip",compression_opts=9,shuffle=True)
        ds_epix = hits.create_dataset('epix',(0,),dtype=dt_epix,maxshape=(None,),chunks=(500,),compression="gzip",compression_opts=9,shuffle=True)
        # Create a data-set for images metadata:
        dt_imgs=np.dtype([('img',np.dtype('S100')),('ncol',np.dtype('int32')),('nrow',np.dtype('int32')),('ccdncol',np.dtype('int32')),('ccdnrow',np.dtype('int32')),('start',np.dtype('S100')),('end',np.dtype('S100')),('rotime',np.dtype('float32')),('imgid',np.dtype('int32'))])
        ds_imgs=f.create_dataset('imgs_metadata',(0,),dtype=dt_imgs,maxshape=(None,),chunks=(500,),compression="gzip",compression_opts=9,shuffle=True)
    return 0

def fillHDF5(filename,events):
    ne=events[0].size
    with h5py.File(filename, "r+") as f:
        ds_param=f['hits/parameters']
        ds_xpix=f['hits/xpix']
        ds_ypix=f['hits/ypix']
        ds_pix=f['hits/pix']
        ds_epix=f['hits/epix']
        # Resize HDF5
        ds_param.resize((ds_param.shape[0]+ne,)) 
        ds_xpix.resize((ds_xpix.shape[0]+ne,)) 
        ds_ypix.resize((ds_ypix.shape[0]+ne,)) 
        ds_pix.resize((ds_pix.shape[0]+ne,)) 
        ds_epix.resize((ds_epix.shape[0]+ne,)) 
        # Fill HDF5
        ds_param['imgid',-ne:]=events[0]
        ds_param['ohdu',-ne:]=events[1]
        ds_param['flag',-ne:]=events[2]
        ds_param['xmin',-ne:]=events[3]
        ds_param['xmax',-ne:]=events[4]
        ds_param['ymin',-ne:]=events[5]
        ds_param['ymax',-ne:]=events[6]
        ds_param['q',-ne:]=events[7]
        ds_param['qe',-ne:]=events[8]
        ds_param['qq',-ne:]=events[9]
        ds_param['npix',-ne:]=events[10]
        ds_param['xbary',-ne:]=events[11]
        ds_param['ybary',-ne:]=events[12]
        ds_param['xvar',-ne:]=events[13]
        ds_param['yvar',-ne:]=events[14]
        ds_xpix[-ne:]=events[15]
        ds_ypix[-ne:]=events[16]
        ds_pix[-ne:]=events[17]
        ds_epix[-ne:]=events[18]
    return 0

def fillHDF5_imgsmtdt(filename,imgsmtdt):
    with h5py.File(filename, "r+") as f:
        ds_imgs=f['imgs_metadata']
        for img in imgsmtdt:
            ds_imgs.resize((ds_imgs.shape[0]+1,))
            ds_imgs['img',-1]=img[0]
            ds_imgs['ncol',-1]=img[1]
            ds_imgs['nrow',-1]=img[2]
            ds_imgs['ccdncol',-1]=img[3]
            ds_imgs['ccdnrow',-1]=img[4]
            ds_imgs['start',-1]=img[5]
            ds_imgs['end',-1]=img[6]
            ds_imgs['rotime',-1]=img[7]
            ds_imgs['imgid',-1]=img[8]
    return 0

def extract_hdu(imgid,hdul,hdu,T1,T2,gain,satvalue):
    print("\nDoing HDU #",hdu)        
    image = hdul[hdu].data
    # Check for NaN pixels
    print(np.sum(np.isnan(image))," NaN pixels")
    # Extract
    label_im, nb_labels = ndimage.label(image>=T2,structure=[[1,1,1],[1,1,1],[1,1,1]])
    labels=np.arange(1,nb_labels+1,dtype=label_im.dtype)
    labelmax=ndimage.maximum(image,labels=label_im,index=labels) # Maximum of each region.
    labelsT1=labels[labelmax>=T1] # Regions that have at least one pixel with value higher than T1.
    nevents=labelsT1.size
    print(nevents," events")
    # Analysis of each region (event):
    iid=np.zeros(nevents)
    ohdu=np.zeros(nevents)
    flag=np.zeros(nevents)
    xMin=np.zeros(nevents)
    xMax=np.zeros(nevents)
    yMin=np.zeros(nevents)
    yMax=np.zeros(nevents)
    e=np.zeros(nevents)  # (adu) event charge
    ee=np.zeros(nevents) # (electrons) event charge
    eq=np.zeros(nevents) # (electrons) quantized event charge
    n=np.zeros(nevents)
    xBary=np.zeros(nevents)
    yBary=np.zeros(nevents)
    xVar=np.zeros(nevents)
    yVar=np.zeros(nevents)
    xPix=[]
    yPix=[]
    ePix=[]
    eePix=[]
    for lid in range(0,labelsT1.size): # Loop over each extracted event
        xy=np.where(label_im==(lid+1))
        pixels=image[xy]
        #Fill the numpy arrays:
        iid[lid]=imgid
        ohdu[lid]=hdu
        flag[lid]=doflag(pixels,satvalue)
        xMin[lid]=np.min(xy[1])
        xMax[lid]=np.max(xy[1])
        yMin[lid]=np.min(xy[0])
        yMax[lid]=np.max(xy[0])
        e[lid]=np.sum(pixels) # Event charge in ADU
        ee[lid]=np.sum(pixels/gain) # Event charge in electrons
        eq[lid]=np.sum(np.around(pixels/gain)) # Quantized event charge 
        n[lid]=pixels.size
        xBary[lid]=np.sum(xy[1]*pixels)/e[lid]
        yBary[lid]=np.sum(xy[0]*pixels)/e[lid]
        xVar[lid]=np.sum((xy[1]-xBary[lid])*(xy[1]-xBary[lid])*pixels)/e[lid]
        yVar[lid]=np.sum((xy[0]-yBary[lid])*(xy[0]-yBary[lid])*pixels)/e[lid]
        xPix.append(xy[1])
        yPix.append(xy[0])
        ePix.append(pixels)
        eePix.append(pixels/gain)
    return [iid,ohdu,flag,xMin,xMax,yMin,yMax,e,ee,eq,n,xBary,yBary,xVar,yVar,xPix,yPix,ePix,eePix]

def doflag(pixels,satvalue):
    flag=int('00000000',2)
    # Check NaN values
    if np.sum(np.isnan(pixels))>0:
        flag = flag | int('00000001',2)
    # Check saturated pixels
    if pixels[pixels>satvalue].size > 0:
        flag = flag | int('00000010',2)
    return flag

def getThr(g,n,i):
    # Skipper-CCD:
    if math.isnan(g[i]):
        # Regular-CCD:
        return n[i]*4.0,n[i]*3.0
    else:
        T1=g[i]/2.0
        T2=g[i]/2.0
        return T1,T2

def initOutFile(outfile):
    if os.path.exists(outfile):
        print("The output file already exist.")
        if query_yes_no("Do you want to overwrite the file?"):
            initHDF5(outfile)
        else:
            print("Bye.")
            exit()
    else:
        initHDF5(outfile)
    return 0

def get_images(folder):
    images=[]
    for filename in os.listdir(folder):
        if filename.endswith(".fits") or filename.endswith(".fits.bz2"):
            images.append(folder+"/"+filename)
    # Extract adquisition time from image header
    imgsmtdt=[] # images metadata.
    for i in images:
        hdul = fits.open(i)
        ncol=hdul[0].header['NCOL']
        nrow=hdul[0].header['NROW']
        ccdncol=hdul[0].header['CCDNCOL']
        ccdnrow=hdul[0].header['CCDNROW']
        dateStart=datetime.strptime(hdul[0].header['DATESTART'],"%Y-%m-%dT%H:%M:%S")
        dateEnd=datetime.strptime(hdul[0].header['DATEEND'],"%Y-%m-%dT%H:%M:%S")
        rotime=(dateEnd-dateStart).total_seconds()
        imgsmtdt.append([i,ncol,nrow,ccdncol,ccdnrow,dateStart,dateEnd,rotime])
        hdul.close()
    imgsmtdt.sort(key=lambda x : x[5]) # Sort the images (older first)
    for x in range(0,len(imgsmtdt)): # Append image ID
        imgsmtdt[x].append(x)
    return imgsmtdt

def query_yes_no(question, default="no"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")
    return 0
