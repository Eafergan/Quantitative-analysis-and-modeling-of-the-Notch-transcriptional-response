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




path=('D:\\AOS\\2nd_AOS_transient_100RBPJ_30Cer\\2nd_AOS_transient_100RBPJ_30Cer_A\\input\\')

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
# for x in onlyfiles:
#     rawfile = bioformats.ImageReader(join(path,x))
#     print (x)
#     for t in range(NumOfT):
#         if t%6==0:
#             print(t)
#         tempimg=rawfile.read(c=0,z=0,t=t,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
#         timeToFile=int((t-(t%2))/2)      
#         if t%2==1:
#             bioformats.write_image(path[0:-6] + 'output\\channel2\\'+x[:-8]+'.tiff', tempimg, PixelT, c=0, z=0, t=timeToFile, size_c=1, size_z=1, size_t=NumOfTReal, channel_names=None)
#         else:
#              bioformats.write_image(path[0:-6] + 'output\\channel1\\'+x[:-8]+'.tiff', tempimg, PixelT, c=0, z=0, t=timeToFile, size_c=1, size_z=1, size_t=NumOfTReal, channel_names=None)          
#         del tempimg, 
#     rawfile.close()
#     del rawfile 
# del x


# onlyfilesFullPath=['D:/AOS/2nd_AOS_transient_100RBPJ_30Cer/2nd_AOS_transient_100RBPJ_30Cer_E/output/channel1/'+string[:-8] + '.tiff' for string in onlyfiles ]
# filesCH1 = " ".join(f'"{f}"' for f in onlyfilesFullPath)

# onlyfilesFullPath=['D:/AOS/2nd_AOS_transient_100RBPJ_30Cer/2nd_AOS_transient_100RBPJ_30Cer_E/output/channel2/'+string[:-8] + '.tiff' for string in onlyfiles ]
# filesCH2 = " ".join(f'"{f}"' for f in onlyfilesFullPath)

# os.chdir("C:\\Program Files\\ilastik-1.4.0rc8-gpu")
# commandCh1 = 'ilastik.exe --headless --project="D:/AOS/2nd_AOS_transient_100RBPJ_30Cer/Ch1_Pxl_n_Obj.ilp" --export_source="Object Predictions"  'f' {filesCH1}'
# os.system(commandCh1)

# commandCh2 = 'ilastik.exe --headless --project="D:/AOS/2nd_AOS_transient_100RBPJ_30Cer/Ch2_Notch4bg_100RBPJ_30Cer.ilp" --export_source="Simple Segmentation"  'f' {filesCH2}'
# os.system(commandCh2)




resultsCh2=np.zeros((NumOfTReal,len(onlyfiles) )) 
bgCh2=np.zeros((NumOfTReal,len(onlyfiles) )) 
MaskCh2=np.zeros((NumOfTReal,len(onlyfiles) )) 


for x in range(len(onlyfiles)):
    nameoffile=  onlyfiles[x]  
    print(nameoffile[:-8]) 
    Ilastik = h5py.File(path[0:-6] + 'output\\channel1\\Ilastik\\'+nameoffile[:-8]+'.h5')
    Mask = Ilastik['exported_data']
    IlastikNotchBG = h5py.File(path[0:-6] + 'output\\channel2\\Ilastik\\'+nameoffile[:-8]+'.h5')
    MaskNotchBG = IlastikNotchBG['exported_data']
    
    for t in range(NumOfTReal):
        tempmask=Mask[t,:,:,0]
        tempbgmaskNotch=MaskNotchBG[t,:,:,0]
        rawfile2 = bioformats.ImageReader( path[0:-6] + 'output\\channel2\\'+nameoffile[:-8]+'.tiff')
        tempmeasure2=rawfile2.read(c=0,z=0,t=t,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)

        tempmask=tempmask*(tempmask==2)/2
        tempbgmask=tempbgmaskNotch*(tempbgmaskNotch==1)*1
        tempNotchMask=tempbgmaskNotch*(tempbgmaskNotch==2)/2

        measureMasked2=tempmeasure2*tempmask
        bgMasked2=tempbgmask*tempmeasure2


        SigCh2=measureMasked2[(measureMasked2<400) & (measureMasked2>0)] 
        BgCh2=bgMasked2[(bgMasked2<400) & (bgMasked2>0)]

        averageSig2=np.sum(np.float64(SigCh2))/len(SigCh2)
        averageBg2=np.sum(np.float64(BgCh2))/len(BgCh2)

        resultsCh2[t,x]=averageSig2-averageBg2

        bgCh2[t,x]=averageBg2
        MaskCh2[t,x]=np.sum(np.float64(tempmask))



        bioformats.write_image(path[0:-6] + 'output\\controls\\'+nameoffile[:-8]+'__measureCh2.tif', measureMasked2,  PixelT, c=0, z=0, t=t, size_c=1, size_z=1, size_t=NumOfTReal, channel_names=None)
        bioformats.write_image(path[0:-6] + 'output\\controls\\'+nameoffile[:-8]+'__bgCh2.tif', bgMasked2,  PixelT, c=0, z=0, t=t, size_c=1, size_z=1, size_t=NumOfTReal, channel_names=None)
        #bioformats.write_image(path[0:-6] + 'output\\controls\\'+nameoffile[:-8]+'__aLLnOTCH.tif', measureAllNotch2,  PixelT, c=0, z=0, t=t, size_c=1, size_z=1, size_t=NumOfTReal, channel_names=None)

        if t%6==0:
            print(t)
        del  tempbgmaskNotch, tempbgmask,  rawfile2, tempmask, tempmeasure2,   measureMasked2,   averageSig2, averageBg2,  bgMasked2,  SigCh2,  BgCh2
    del nameoffile, Mask, Ilastik


NamesVector = [elem[elem.find('__')+2:elem.find('.ome')] for elem in onlyfiles]




