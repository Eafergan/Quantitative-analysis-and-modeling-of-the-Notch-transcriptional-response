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



path=('D:\\Zeiss\\TSA2uM\\input\\')

from os.path import isfile, join
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

#assuming there is no more than 100 hubs , 0:size 1:JFsignal 2:CerSignal 3: num of hubs 4: average nuclear JF sig 5: average nuclear Cer Sig
Results=np.zeros((len(onlyfiles),100,6)) 
Scatter=np.zeros((2000,3)) 
Norm_Scatter=np.zeros((2000,3)) 

counter=0 #to get a single vector of all hub

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
    rawfileMaskJustNucleoli = bioformats.ImageReader(join(path[:-6]+'output\\channel2_Cer\\Nucleoli\\' + nameOfFile+ 'f'))



    print (nameOfFile)
    Mask=np.zeros((344,344,NumOfZ))
    JF=np.zeros((344,344,NumOfZ))
    Cer=np.zeros((344,344,NumOfZ))
    MaskJustNucleus=np.zeros((344,344,NumOfZ))
    MaskJustNucleoli=np.zeros((344,344,NumOfZ))

    for z in range(NumOfZ):  #Ilastik marked the Z slices as T slices, this part loads the images along all z slices
        Mask[:,:,z]=rawfileMask.read(c=0,z=0,t=z,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        Cer[:,:,z]=rawfileCer.read(c=0,z=z,t=0,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        JF[:,:,z]=rawfilecleanJF.read(c=0,z=z,t=0,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        MaskJustNucleus[:,:,z]=rawfileMaskJustNucleus.read(c=0,z=z,t=0,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        MaskJustNucleoli[:,:,z]=rawfileMaskJustNucleoli.read(c=0,z=0,t=z,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)

    rawfileMask.close()
    rawfilecleanJF.close()
    rawfileCer.close()
    rawfileMaskJustNucleus.close()
    rawfileMaskJustNucleoli.close()

    del z
       

    Mask=Mask*(Mask==2)/2 
    labeled_image, num_features = label(Mask)
    
    MaskJustNucleus_woHubs=MaskJustNucleus-Mask
    MaskJustNucleus_woHubs=MaskJustNucleus_woHubs*(MaskJustNucleus_woHubs==1)/1

    MaskJustNucleoli=MaskJustNucleoli*(MaskJustNucleoli==3)/3
    MaskJustNucleus_woNucleoli=MaskJustNucleus_woHubs-MaskJustNucleoli
    MaskJustNucleus_woNucleoli=MaskJustNucleus_woNucleoli*(MaskJustNucleus_woNucleoli==1)/1
    
    CleanCer=Cer*MaskJustNucleus_woNucleoli
    CleanJF=JF*MaskJustNucleus_woNucleoli
    
    CerInNucleoli=MaskJustNucleoli*Cer
    JF_in_nucleoli=MaskJustNucleoli*JF
    
    Results[x,0,4]=np.sum(CleanJF)/np.sum(MaskJustNucleus_woNucleoli) #mean signal without hubs and nocleoli
    Results[x,0,5]=np.sum(CleanCer)/np.sum(MaskJustNucleus_woNucleoli)
    Results[x,1,4]=np.sum(JF_in_nucleoli)/np.sum(MaskJustNucleoli) #mean signal in nucleoli
    Results[x,1,5]=np.sum(CerInNucleoli)/np.sum(MaskJustNucleoli)
    Results[x,3,4]=np.sum(MaskJustNucleoli)
    Results[x,3,5]=np.sum(MaskJustNucleus_woNucleoli)
    
    AllHubsJF=Mask*JF
    AllHubsCer=Mask*Cer
    
    EntireCellSignalJF=MaskJustNucleus*JF
    EntireCellSignalCer=MaskJustNucleus*Cer

    
    Results[x,2,4]=np.sum(AllHubsJF)/np.sum(EntireCellSignalJF) #mean signal in nucleoli
    Results[x,2,5]=np.sum(AllHubsCer)/np.sum(EntireCellSignalCer)

    del AllHubsJF,AllHubsCer, EntireCellSignalJF, EntireCellSignalCer
    
    
    
    for i in range(1,num_features+1):
    #    print(i)
        Results[x,i,3]=i
        single_hub=labeled_image*(labeled_image==i)/i
        iCer=Cer*single_hub
        iJF=JF*single_hub
        Results[x,i,0]=np.sum(single_hub) #volume
        Results[x,i,1]=np.sum(iJF)/np.sum(single_hub) #mean signal
        Results[x,i,2]=np.sum(iCer)/np.sum(single_hub)
        Scatter[counter,0]=Results[x,i,0]
        Scatter[counter,1]=Results[x,i,1]
        Scatter[counter,2]=Results[x,i,2]
        Norm_Scatter[counter,0]=Results[x,i,0]
        Norm_Scatter[counter,1]=Results[x,i,1]/Results[x,0,4]
        Norm_Scatter[counter,2]=Results[x,i,2]/Results[x,0,5]
        counter+=1
        
    
    
    for z in range(NumOfZ):         #Ilastik marked the Z slices as T slices, this part loads the images along all z slices
        bioformats.write_image(path[0:-6] + 'output\\control\\'+nameOfFile, JF[:,:,z],  PixelT, c=0, z=z, t=0, size_c=2, size_z=NumOfZ, size_t=1, channel_names=None)
        bioformats.write_image(path[0:-6] + 'output\\control\\'+nameOfFile, Mask[:,:,z],  PixelT, c=1, z=z, t=0, size_c=2, size_z=NumOfZ, size_t=1, channel_names=None)

    
    del rawfileMask ,rawfilecleanJF,rawfileMaskJustNucleus, rawfileCer, xml_string, ome, iome, NumOfZ, PixelT, nameOfFile
    del  JF, Cer, MaskJustNucleus, labeled_image, num_features, Mask, 

#plt.imshow(AllHubsJF[:,:,15])
#z=15
#bioformats.write_image(path[0:-6] + 'output\\controlll\\'+nameOfFile, MaskJustNucleoli[:,:,z],  PixelT, c=1, z=z, t=0, size_c=2, size_z=NumOfZ, size_t=1, channel_names=None)

#np.min(JF_in_nucleoli)






    
    
    