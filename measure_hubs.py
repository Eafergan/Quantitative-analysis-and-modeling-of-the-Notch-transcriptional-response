%reset -f

import bioformats
import javabridge
javabridge.start_vm(class_path=bioformats.JARS)
import os
from os import listdir
from os.path import isfile, join
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from scipy.ndimage import label
from scipy.ndimage import binary_fill_holes


myloglevel="ERROR"  # user string argument for logLevel.

javabridge.start_vm(class_path=bioformats.JARS)

rootLoggerName = javabridge.get_static_field("org/slf4j/Logger","ROOT_LOGGER_NAME", "Ljava/lang/String;")
rootLogger = javabridge.static_call("org/slf4j/LoggerFactory","getLogger", "(Ljava/lang/String;)Lorg/slf4j/Logger;", rootLoggerName)
logLevel = javabridge.get_static_field("ch/qos/logback/classic/Level",myloglevel, "Lch/qos/logback/classic/Level;")
javabridge.call(rootLogger, "setLevel", "(Lch/qos/logback/classic/Level;)V", logLevel)



path=('D:\\Zeiss\\measure_RBPJ\\DAPT\\input\\')

from os.path import isfile, join
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

#assuming there is no more than 100 hubs , 0:size 1:JFsignal 2:CerSignal 3: num of hubs 4: average nuclear JF sig 5: average nuclear Cer Sig
Results=np.zeros((len(onlyfiles),100,6)) 

for x in range(len(onlyfiles)):
    nameOfFile=onlyfiles[x]
    xml_string = bioformats.get_omexml_metadata(path+'\\'+nameOfFile)
    ome = bioformats.OMEXML(xml_string) # be sure everything is ascii
    iome = ome.image(0) # e.g. first image
    NumOfZ=iome.Pixels.get_SizeZ()
    PixelT=iome.Pixels.get_PixelType()
    rawfilecleanJF = bioformats.ImageReader(join(path[:-6]+'output\\cleaned\\' + nameOfFile))
    rawfileCer = bioformats.ImageReader(join(path[:-6]+'output\\channel2_Cer\\' + nameOfFile+ 'f'))
    rawfileMask = bioformats.ImageReader(join(path[:-6]+'output\\cleaned\\Ilastik\\' + nameOfFile+ 'f'))
    rawfileMaskJustNucleus = bioformats.ImageReader(join(path[:-6]+'output\\cleanedMask\\' + nameOfFile))

    print (nameOfFile)
    Mask=np.zeros((344,344,NumOfZ))
    JF=np.zeros((344,344,NumOfZ))
    Cer=np.zeros((344,344,NumOfZ))
    MaskJustNucleus=np.zeros((344,344,NumOfZ))

    for z in range(NumOfZ):  #Ilastik marked the Z slices as T slices, this part loads the images along all z slices
        Mask[:,:,z]=rawfileMask.read(c=0,z=0,t=z,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        Cer[:,:,z]=rawfileCer.read(c=0,z=z,t=0,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        JF[:,:,z]=rawfilecleanJF.read(c=0,z=z,t=0,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        MaskJustNucleus[:,:,z]=rawfileMaskJustNucleus.read(c=0,z=z,t=0,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
   
    rawfileMask.close()
    rawfilecleanJF.close()
    rawfileCer.close()
    rawfileMaskJustNucleus.close()
    del z
       

    Mask=Mask*(Mask==2)/2 
    labeled_image, num_features = label(Mask)
    
    MaskJustNucleus_woHubs=MaskJustNucleus-Mask
    MaskJustNucleus_woHubs=MaskJustNucleus_woHubs*(MaskJustNucleus_woHubs==1)/1

    
    CleanCer=Cer*MaskJustNucleus_woHubs
    CleanJF=JF*MaskJustNucleus_woHubs
    Results[x,0,4]=np.sum(CleanJF)/np.sum(MaskJustNucleus_woHubs)
    Results[x,0,5]=np.sum(CleanCer)/np.sum(MaskJustNucleus_woHubs)

    
    for i in range(1,num_features+1):
    #    print(i)
        Results[x,i,3]=i
        single_hub=labeled_image*(labeled_image==i)/i
        iCer=Cer*single_hub
        iJF=JF*single_hub
        Results[x,i,0]=np.sum(single_hub)
        Results[x,i,1]=np.sum(iJF)/np.sum(single_hub)
        Results[x,i,2]=np.sum(iCer)/np.sum(single_hub)
    
    
 #   for z in range(NumOfZ):  #Ilastik marked the Z slices as T slices, this part loads the images along all z slices
 #       bioformats.write_image(path[0:-6] + 'output\\control\\'+nameOfFile, JF[:,:,z],  PixelT, c=0, z=z, t=0, size_c=2, size_z=NumOfZ, size_t=1, channel_names=None)
   #     bioformats.write_image(path[0:-6] + 'output\\control\\'+nameOfFile, Mask[:,:,z],  PixelT, c=1, z=z, t=0, size_c=2, size_z=NumOfZ, size_t=1, channel_names=None)

    
    del rawfileMask ,rawfilecleanJF,rawfileMaskJustNucleus, rawfileCer, xml_string, ome, iome, NumOfZ, PixelT, nameOfFile
    del  JF, Cer, MaskJustNucleus, labeled_image, num_features, Mask, MaskJustNucleus_woHubs

plt.imshow(CleanJF[:,:,16])


#np.max(MaskJustNucleus)





    
    
    