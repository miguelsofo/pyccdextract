import numpy as np
import matplotlib.pyplot as plt 
from scipy import stats
from scipy import ndimage
import ccd


rn=2
T2=1.5*rn
image=ccd.simNoiseImg(rn,10,10)
ccd.plot_fits(image,T2,10)

label_im, nb_labels = ndimage.label(image>=T2,structure=[[1,1,1],[1,1,1],[1,1,1]])
print(label_im)
plt.figure()
plt.imshow(label_im)
plt.show()

labels=np.arange(1,nb_labels+1,dtype=label_im.dtype)
labelmax=ndimage.maximum(image,labels=label_im,index=labels) # Maximum of each region.
print(np.where(label_im==labels))

labelsT1=labels[labelmax>=T1] # Regions that have at least one pixel with value higher than T1.
nevents=labelsT1.size
