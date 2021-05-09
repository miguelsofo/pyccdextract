from skimage.measure import label, regionprops, regionprops_table
from skimage.filters import apply_hysteresis_threshold
from skimage.segmentation import clear_border

#from scipy import ndimage
from astropy.io import fits
from datetime import datetime
import matplotlib.pyplot as plt 
import numpy as np
import math
import sys
import h5py
import os.path
import pandas as pd

def runExtract(config):
    # -----------------------------------------------------------     
    try:
        outfile=config["out_file"] # Output file.
    except:
        print('Please provide an output file name. Bye.')
        exit()
    initOutFile(outfile) # Check and init the output HDF5 file.
    # -----------------------------------------------------------         
    try: 
        extinput=config["input"]
    except:
        print('Please provide an input. Bye')
        exit()
    idf=read_folder(config) # Generates images metadata.
    idf.to_csv(outfile.rsplit('/',1)[0]+'/imgs_metadata.csv')
    # -----------------------------------------------------------     
    # Start extraction process: 
    # Loop over the images: 
    for i in idf.index:
        imgname=idf['img'][i]
        imgid=idf['imgid'][i] # image ID from image metadata.
        print("\n===============================================================")
        print(imgname)
        hdul = fits.open(imgname)
        for hdu in range(0,len(hdul)):
            T1=idf['T1_'+str(hdu)][i]
            T2=idf['T2_'+str(hdu)][i]
            gain=idf['gain_'+str(hdu)][i]
            satvalue=idf['satValue'][i]
            events=extract_hdu(imgid,hdul,hdu,T1,T2,gain,satvalue)
            fillHDF5(outfile,events)
    fillHDF5_imgsmtdt(outfile,idf) # Fill HDF5 with images metadata.
    return 0

def get_satvalue(imgname):
    hdul = fits.open(imgname) # read the image
    hdumax=np.zeros(len(hdul))
    for hdu in range(0,len(hdul)): # loop over image hdus
        tmpmax=np.max(hdul[hdu].data)
        if tmpmax>hdumax[hdu]:
            hdumax[hdu]=tmpmax
    hdumax=hdumax*0.95
    return np.min(hdumax)

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
        #dt_imgs=np.dtype([('img',np.dtype('S100')),('ncol',np.dtype('int32')),('nrow',np.dtype('int32')),('ccdncol',np.dtype('int32')),('ccdnrow',np.dtype('int32')),('start',np.dtype('S100')),('end',np.dtype('S100')),('rotime',np.dtype('float32')),('imgid',np.dtype('int32'))])
        #ds_imgs=f.create_dataset('imgs_metadata',(0,),dtype=dt_imgs,maxshape=(None,),chunks=(500,),compression="gzip",compression_opts=9,shuffle=True)
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
    hdf=pd.HDFStore(filename,'r+')
    hdf.put('/imgs_metadata',imgsmtdt)
    return 0

def doflag(pixels,satvalue):
    flag=int('00000000',2)
    # Check NaN values
    if np.sum(np.isnan(pixels))>0:
        flag = flag | int('00000001',2)
    # Check saturated pixels
    if pixels[pixels>satvalue].size > 0:
        flag = flag | int('00000010',2)
    return flag

def getThr(g,n):
    # Skipper-CCD:
    if math.isnan(g):
        # Regular-CCD:
        return n*4.0,n*3.0
    else:
        # Skipper-CCD
        return g/2.0,g/2.0

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

def check_imgs(images):
    # This routine check if all images have the same number of hdu.
    nhdu=np.zeros(len(images))
    for i in range(0,len(images)):
        nhdu[i]=len(fits.open(images[i]))
    res=np.all(nhdu==nhdu[0])
    if res:
        print("All images have %d HDUs."%nhdu[0])
    else:
        print("The images does not have the same number of HDUs. Not allowed")
        exit()
    return int(nhdu[0])

def read_folder(config):
    folder=config['input']
    if not os.path.exists(folder):
        print("The input folder does not exist. Bye...")
        exit()
    # -----------------------------------------------------------    
    # List of images:
    images=[]
    for filename in os.listdir(folder):
        if filename.endswith(".fits") or filename.endswith(".fits.bz2"):
            images.append(folder+"/"+filename)
    nhdu=check_imgs(images) # Check if all images have the same number of HDUs
    # -----------------------------------------------------------    
    try:
        n=config["noise"] # Noise of each CCD amplifier (extension). 
    except:
        print("Please provide a vector with the noise of each HDU.")
        exit()
    if(len(n)!=nhdu):
        print("The vector length with the noises does not match the number of HDUs of the images.")
        exit()
    # -----------------------------------------------------------     
    try: 
        g=config["gain"] # (ADU/e) Gain of each CCD amplifier (extension).
        print("Extraction of sub-electron noise images.")
    except:
        print("GAIN values not provied -> Extraction of regular-CCD images.")
        g=np.ones(len(n))*np.nan
    if(len(g)!=nhdu):
        print("The vector length with the gains does not match the number of HDUs of the images.")
        exit()
    # -----------------------------------------------------------     
    hdr=[]
    try:
        hdr=config['hdr']
    except:
        print("No additional header words provided.")
    # -----------------------------------------------------------     
    imgs_mtdt=get_imgsdf(images,hdr)  
    # -----------------------------------------------------------     
    # Append additional data necessary during extraction.
    for i in range(nhdu):
        imgs_mtdt["std_"+str(i)]=np.ones(len(imgs_mtdt['img']))*n[i]
        imgs_mtdt["gain_"+str(i)]=np.ones(len(imgs_mtdt['img']))*g[i]
        T1,T2=getThr(g[i],n[i])
        imgs_mtdt["T1_"+str(i)]=np.ones(len(imgs_mtdt['img']))*T1
        imgs_mtdt["T2_"+str(i)]=np.ones(len(imgs_mtdt['img']))*T2
    return imgs_mtdt

def get_imgsdf(images,hdr=[]):
    # images: list of images to process. 
    # hdr: list of additional header words that you want to include in the images metadata.
    # Init dataframe:
    column_names=['img','imgid','ncol','nrow','ccdncol','ccdnrow','dateStart','dateEnd','roTime','pixrot','satValue']
    if len(hdr)>0:
        column_names.append(hdr)
    df=pd.DataFrame(columns=column_names)
    # Extract adquisition time from image header
    for i in range(len(images)):
        imgname=images[i]
        hdul = fits.open(imgname)
        df.loc[i,'img']=imgname
        df.loc[i,'imgid']=i
        df.loc[i,'ncol']=float(hdul[0].header['NCOL'])
        df.loc[i,'nrow']=float(hdul[0].header['NROW'])
        df.loc[i,'ccdncol']=float(hdul[0].header['CCDNCOL'])
        df.loc[i,'ccdnrow']=float(hdul[0].header['CCDNROW'])
        df.loc[i,'dateStart']=datetime.strptime(hdul[0].header['DATESTART'],"%Y-%m-%dT%H:%M:%S")
        df.loc[i,'dateEnd']=datetime.strptime(hdul[0].header['DATEEND'],"%Y-%m-%dT%H:%M:%S")
        df.loc[i,'roTime']=(df.loc[i,'dateEnd']-df.loc[i,'dateStart']).total_seconds()
        df.loc[i,'pixrot']=(df.loc[i,'nrow']*1e3)/(df.loc[i,'ncol']*df.loc[i,'nrow'])
        df.loc[i,'satValue']=get_satvalue(df.loc[i,'img'])
        # Add additional headers:
        for h in hdr:
            df.loc[i,h]=hdul[0].header[h]
        hdul.close()
    df.info()
    return df

def extract_hdu(imgid,hdul,hdu,T1,T2,gain,satvalue):
    print("\nDoing HDU #",hdu)        
    image = hdul[hdu].data
    # Check for NaN pixels
    print(np.sum(np.isnan(image))," NaN pixels")
    # get the mask via hysteresis threshold and remove events in borders
    mask = apply_hysteresis_threshold(image,float(T2),float(T1))
    #clear_border(mask, in_place=True)
    # get the labels
    label_im, nevents = label(mask, connectivity=2, return_num=True)
    print('> Found {} different regions'.format(nevents))
    # Analysis of each region (event):
    props = regionprops_table(label_im, intensity_image=image, properties=('area', 'coords', 'intensity_image', 'label', 'weighted_centroid', 'weighted_moments_central', 'bbox'))
    # the important vectors
    iid = props['label']
    ohdu = np.ones(nevents)*hdu
    xMin = props['bbox-1']
    xMax = props['bbox-3']-1
    yMin = props['bbox-0']
    yMax = props['bbox-2']-1
    e = props['weighted_moments_central-0-0']
    ee = e/gain
    n = props['area']
    xBary = props['weighted_centroid-1']
    yBary = props['weighted_centroid-0']
    xVar = props['weighted_moments_central-0-2']/e
    yVar = props['weighted_moments_central-2-0']/e
    yPix, xPix = [ [ props['coords'][i][:,j] for i in range(nevents)] for j in [0,1]]
    ePix = [ image[(yPix[i], xPix[i])] for i in range(nevents) ]
    eePix = [ ePix[i]/gain for i in range(nevents) ]
    flag = np.array([doflag(ePix[i],satvalue) for i in range(nevents)])
    eq = np.array( [np.around(ePix[i]/gain).sum() for i in range(nevents)])
    # all done
    return [iid,ohdu,flag,xMin,xMax,yMin,yMax,e,ee,eq,n,xBary,yBary,xVar,yVar,xPix,yPix,ePix,eePix]

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
