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


x=onlyfiles[1]

for x in onlyfiles:
    xml_string = bioformats.get_omexml_metadata(path+'\\'+x)
    ome = bioformats.OMEXML(xml_string) # be sure everything is ascii
    iome = ome.image(0) # e.g. first image
    NumOfZ=iome.Pixels.get_SizeZ()
    PixelT=iome.Pixels.get_PixelType()
    rawfileJF = bioformats.ImageReader(join(path[:-6]+'output\\channel1_JF\\' + x+ 'f'))
    rawfileCer = bioformats.ImageReader(join(path[:-6]+'output\\channel2_Cer\\' + x+ 'f'))
    rawfileMask = bioformats.ImageReader(join(path[:-6]+'output\\channel2_Cer\\Ilastik\\' + x+ 'f'))
    print (x)
    Mask=np.zeros((344,344,NumOfZ))
    JF=np.zeros((344,344,NumOfZ))
    
    for z in range(NumOfZ):  #Ilastik marked the Z slices as T slices, this part loads the images along all z slices
        Mask[:,:,z]=rawfileMask.read(c=0,z=0,t=z,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        JF[:,:,z]=rawfileJF.read(c=0,z=z,t=0,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
    rawfileMask.close()
    rawfileJF.close()
    del z
    
    tempmask=Mask
    tempmask=tempmask*(tempmask==3)/3 #this creates an image of the ilastik nuclues preditction alone
    labeled_image, num_features = label(tempmask)
    
    # Find the sizes of each connected component
    component_sizes = np.bincount(labeled_image.ravel())
    # Exclude the zero component (background) and find the largest non-zero component
    non_zero_component_sizes = component_sizes[1:]
    largest_component = non_zero_component_sizes.argmax() + 1
    nucleus=labeled_image*(labeled_image==largest_component)/largest_component
    del labeled_image, component_sizes, num_features, tempmask, non_zero_component_sizes, largest_component
    
    filled_nucleus=np.zeros((344,344,NumOfZ))
    for z in range(NumOfZ):
        filled_nucleus[:,:,z] = binary_fill_holes(nucleus[:,:,z])
    Clean_JF=JF*filled_nucleus
    
    for z in range(NumOfZ):  #Ilastik marked the Z slices as T slices, this part loads the images along all z slices
        bioformats.write_image(path[0:-6] + 'output\\cleaned\\'+x, Clean_JF[:,:,z],  PixelT, c=0, z=z, t=0, size_c=1, size_z=NumOfZ, size_t=1, channel_names=None)
        bioformats.write_image(path[0:-6] + 'output\\cleanedMask\\'+x, filled_nucleus[:,:,z],  PixelT, c=0, z=z, t=0, size_c=1, size_z=NumOfZ, size_t=1, channel_names=None)

    
    del rawfileJF,rawfileMask , xml_string, ome, iome, NumOfZ, PixelT, filled_nucleus,Clean_JF

 







    
    
    