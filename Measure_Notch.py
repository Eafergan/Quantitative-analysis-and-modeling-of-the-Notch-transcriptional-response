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
import cv2 
import scipy.ndimage as ndimage




myloglevel="ERROR"  # user string argument for logLevel.

javabridge.start_vm(class_path=bioformats.JARS)

rootLoggerName = javabridge.get_static_field("org/slf4j/Logger","ROOT_LOGGER_NAME", "Ljava/lang/String;")
rootLogger = javabridge.static_call("org/slf4j/LoggerFactory","getLogger", "(Ljava/lang/String;)Lorg/slf4j/Logger;", rootLoggerName)
logLevel = javabridge.get_static_field("ch/qos/logback/classic/Level",myloglevel, "Lch/qos/logback/classic/Level;")
javabridge.call(rootLogger, "setLevel", "(Lch/qos/logback/classic/Level;)V", logLevel)



path=('D:\\Zeiss\\measure_Notchs\\Notch+\\input\\')

from os.path import isfile, join
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

#assuming there is no more than 100 hubs , 0:size 1:JFsignal 2:CerSignal 3: num of hubs 4: average nuclear JF sig 5: average nuclear Cer Sig
Results=np.zeros((len(onlyfiles),3)) 


for x in range(len(onlyfiles)):
    nameOfFile=onlyfiles[x]
    xml_string = bioformats.get_omexml_metadata(path+'\\'+nameOfFile)
    ome = bioformats.OMEXML(xml_string) # be sure everything is ascii
    iome = ome.image(0) # e.g. first image
    NumOfZ=iome.Pixels.get_SizeZ()
    PixelT=iome.Pixels.get_PixelType()
    rawfileJF = bioformats.ImageReader(join(path[:-6]+'output\\channel1_JF\\' + nameOfFile+ 'f'))
    rawfileCer = bioformats.ImageReader(join(path[:-6]+'output\\channel2_Cer\\' + nameOfFile+ 'f'))
    rawfileMask = bioformats.ImageReader(join(path[:-6]+'output\\channel2_Cer\\Ilastik\\' + nameOfFile+ 'f'))


    print (nameOfFile)
    Mask=np.zeros((344,344,NumOfZ))
    JF=np.zeros((344,344,NumOfZ))
    Cer=np.zeros((344,344,NumOfZ))

    for z in range(NumOfZ):  #Ilastik marked the Z slices as T slices, this part loads the images along all z slices
        Mask[:,:,z]=rawfileMask.read(c=0,z=0,t=z,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        Cer[:,:,z]=rawfileCer.read(c=0,z=z,t=0,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        JF[:,:,z]=rawfileJF.read(c=0,z=z,t=0,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)

    rawfileMask.close()
    rawfileJF.close()
    rawfileCer.close()

    
    del z
    
    Mask=Mask*(Mask==2)/2 


    
    
    filled_nucleus=np.zeros((344,344,NumOfZ))
    for z in range(NumOfZ):
        filled_nucleus[:,:,z] = binary_fill_holes(Mask[:,:,z])

    
    labeled_image, num_features = label(filled_nucleus)
    
    # Find the sizes of each connected component
    component_sizes = np.bincount(labeled_image.ravel())
    # Exclude the zero component (background) and find the largest non-zero component
    non_zero_component_sizes = component_sizes[1:]
    largest_component = non_zero_component_sizes.argmax() + 1
    nucleus=labeled_image*(labeled_image==largest_component)/largest_component
    del labeled_image, component_sizes, num_features,  non_zero_component_sizes, largest_component
    
        
    structure_z = np.ones((3, 3, 3), dtype=bool)
   # structure_z[1, 1, 1] = 1  # Only affect the central line along the z-axis

    
    Eroded_nuc = ndimage.binary_erosion(nucleus, structure=structure_z)
    

    CleanCer=Cer*Eroded_nuc
    CleanJF=JF*Eroded_nuc
    
    
    #BackgroundMeasured average between all repeat=1.003    
   # Sig90JF=CleanJF[(CleanJF<np.percentile(CleanJF[CleanJF>0], 90)) & (CleanJF>0)]

    Results[x,0]=np.sum(Eroded_nuc) #Volume of cell for control
    Results[x,1]=(np.sum(CleanJF)/np.sum(Eroded_nuc))
    Results[x,2]=np.sum(CleanCer)/np.sum(Eroded_nuc) #mean signal in nucleoli


    
    for z in range(NumOfZ):         #Ilastik marked the Z slices as T slices, this part loads the images along all z slices
        bioformats.write_image(path[0:-6] + 'output\\control\\'+nameOfFile, CleanJF[:,:,z],  PixelT, c=0, z=z, t=0, size_c=1, size_z=NumOfZ, size_t=1, channel_names=None)

    
    del rawfileMask ,rawfileJF, rawfileCer, xml_string, ome, iome, NumOfZ, PixelT, nameOfFile
    del  JF, Cer,   Mask 

ResultsBoth=np.zeros((200,2)) 


#plt.imshow(JF_BG[:,:,17])
#z=15
#bioformats.write_image(path[0:-6] + 'output\\controlll\\'+nameOfFile, MaskJustNucleoli[:,:,z],  PixelT, c=1, z=z, t=0, size_c=2, size_z=NumOfZ, size_t=1, channel_names=None)

#np.min(JF_in_nucleoli)






    
    
    