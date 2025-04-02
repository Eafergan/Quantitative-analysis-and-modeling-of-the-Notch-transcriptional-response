%reset -f

from matplotlib import pyplot as plt
import numpy as np
import bioformats
from PIL import Image
import cv2
#from skimage.morphology import (erosion, dilation, opening, closing, disk)
import javabridge
javabridge.start_vm(class_path=bioformats.JARS)
from bioformats import ImageReader
import bioformats.omexml as ome
import os
from os import listdir
from os.path import isfile, join
from datetime import (date ,datetime)
import h5py
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import pandas as pd
import subprocess
from pathlib import Path

myloglevel="ERROR"  # user string argument for logLevel.

javabridge.start_vm(class_path=bioformats.JARS)

rootLoggerName = javabridge.get_static_field("org/slf4j/Logger","ROOT_LOGGER_NAME", "Ljava/lang/String;")
rootLogger = javabridge.static_call("org/slf4j/LoggerFactory","getLogger", "(Ljava/lang/String;)Lorg/slf4j/Logger;", rootLoggerName)
logLevel = javabridge.get_static_field("ch/qos/logback/classic/Level",myloglevel, "Lch/qos/logback/classic/Level;")
javabridge.call(rootLogger, "setLevel", "(Lch/qos/logback/classic/Level;)V", logLevel)



def Single_exp_decay(x, C1,C2): 
    return  (1-C2)*np.exp(-x*C1) + C2




path=('D:\\RBPJ\\2nd_RBPJ_DAPTwo_48hrs\\2nd_RBPJ_DAPTwo_48hrs_D\\\input\\')

from os.path import isfile, join
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

nameoffile=onlyfiles[0] # just inorder to have first estimaion of times

xml_string = bioformats.get_omexml_metadata(path+'\\'+nameoffile)
ome = bioformats.OMEXML(xml_string) # be sure everything is ascii
iome = ome.image(0) # e.g. first image
NumOfT=iome.Pixels.get_SizeT()
PixelT=iome.Pixels.get_PixelType()

#this section seperate the raw input file to the 2 channels
NumOfTReal=int(NumOfT/2);
for x in onlyfiles:
    rawfile = bioformats.ImageReader(join(path,x))
    print (x)
    for t in range(NumOfT):
        if t%6==0:
            print(t)
        tempimg=rawfile.read(c=0,z=0,t=t,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        timeToFile=int((t-(t%2))/2)      
        if t%2==1:
            bioformats.write_image(path[0:-6] + 'output\\channel2\\'+x[:-8]+'.tiff', tempimg, PixelT, c=0, z=0, t=timeToFile, size_c=1, size_z=1, size_t=NumOfTReal, channel_names=None)
        else:
             bioformats.write_image(path[0:-6] + 'output\\channel1\\'+x[:-8]+'.tiff', tempimg, PixelT, c=0, z=0, t=timeToFile, size_c=1, size_z=1, size_t=NumOfTReal, channel_names=None)          
        del tempimg, 
    rawfile.close()
    del rawfile 
del x



resultsCh2=np.zeros((NumOfTReal,len(onlyfiles) )) 
bgCh2=np.zeros((NumOfTReal,len(onlyfiles) )) 
MaskCh2=np.zeros((NumOfTReal,len(onlyfiles) )) 


for x in range(len(onlyfiles)):
    nameoffile=  onlyfiles[x]  
    print(nameoffile[:-8]) 
    Ilastik = h5py.File(path[0:-6] + 'output\\channel1\\Ilastik\\'+nameoffile[:-8]+'.h5')
    Mask = Ilastik['exported_data']
 
    for t in range(NumOfTReal):
        tempmask=Mask[t,:,:,0]
        tempbgmask=Mask[t,:,:,0]
        rawfile2 = bioformats.ImageReader( path[0:-6] + 'output\\channel2\\'+nameoffile[:-8]+'.tiff')
        tempmeasure2=rawfile2.read(c=0,z=0,t=t,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)

        tempmask=tempmask*(tempmask==2)/2
        tempbgmask=tempbgmask*(tempbgmask==1)*1
        measureMasked2=tempmeasure2*tempmask

        bgMasked2=tempbgmask*tempmeasure2

        Sig90Ch2=measureMasked2[(measureMasked2<np.percentile(measureMasked2[measureMasked2>0], 90)) & (measureMasked2>0)]
        Bg90Ch2=bgMasked2[(bgMasked2<np.percentile(bgMasked2[bgMasked2>0], 90)) & (bgMasked2>0)]

        averageSig2=np.sum(np.float64(Sig90Ch2))/np.sum(np.float64(tempmask))
        averageBg2=np.sum(np.float64(Bg90Ch2))/np.sum(np.float64(tempbgmask)) ##maybe change size to np.count_nonzero
 
        resultsCh2[t,x]=averageSig2-averageBg2
        bgCh2[t,x]=averageBg2
        MaskCh2[t,x]=np.sum(np.float64(tempmask))



        bioformats.write_image(path[0:-6] + 'output\\controls\\'+nameoffile[:-8]+'__measureCh2.tif', measureMasked2,  PixelT, c=0, z=0, t=t, size_c=1, size_z=1, size_t=NumOfTReal, channel_names=None)
        bioformats.write_image(path[0:-6] + 'output\\controls\\'+nameoffile[:-8]+'__bgCh2.tif', bgMasked2,  PixelT, c=0, z=0, t=t, size_c=1, size_z=1, size_t=NumOfTReal, channel_names=None)
     
  
        if t%6==0:
            print(t)
        del  rawfile2, tempmask, tempmeasure2, measureMasked2, averageSig2, averageBg2, bgMasked2, Sig90Ch2, Bg90Ch2
    del nameoffile, Mask, Ilastik


NamesVector = [elem[elem.find('__')+2:elem.find('.ome')] for elem in onlyfiles]





### so far the results were "resultsCh1", now the for these results to exponential decay.



first_time_points = resultsCh2[0,:]
normalized_matrix = resultsCh2 / first_time_points[np.newaxis, :]

first_time_points_Growth = MaskCh2[0,:]
normalized_matrix_growth =  first_time_points_Growth[np.newaxis, :]/ MaskCh2 





vec = np.linspace(0, (NumOfTReal-1)*20, num=NumOfTReal)

C1Mat=np.zeros((6))
C2Mat=np.zeros((6))
C1MatG=np.zeros((6))
C2MatG=np.zeros((6))


for x in range(6):
    params, covariance = curve_fit(Single_exp_decay, vec, normalized_matrix[:, x],bounds=([0,0],[1,1]))
    y_pred = Single_exp_decay( vec, *params)
    C1Mat[x]=params[0]
    C2Mat[x]=params[1]
    print(x)
    del params, covariance, y_pred 
    
for x in range(6):
    params, covariance = curve_fit(Single_exp_decay, vec, normalized_matrix_growth[:, x],bounds=([0,0],[1,1]))
    y_pred = Single_exp_decay( vec, *params)
    C1MatG[x]=params[0]
    C2MatG[x]=params[1]
    print(x)
    del params, covariance, y_pred 

HalfLife=np.log(2)/C1Mat
HalfLifeGrowth=np.log(2)/C1MatG



for i in range(6):
    y_pred = Single_exp_decay( vec,  C1Mat[i],C2Mat[i])
    plt.plot(vec,y_pred,linestyle='dashed')
    plt.plot(vec,normalized_matrix[ :, i],linestyle='dashed')
    plt.ylim((0,1))
    plt.show()
    
for i in range(6):
    y_pred = Single_exp_decay( vec,  C1MatG[i],C2MatG[i])
    plt.plot(vec,y_pred,linestyle='dashed')
    plt.plot(vec,normalized_matrix_growth[ :, i],linestyle='dashed')
    plt.ylim((0,1))
    plt.show()
    
