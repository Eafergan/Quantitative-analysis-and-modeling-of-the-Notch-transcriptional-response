%reset -f

from matplotlib import pyplot as plt
import numpy as np
import bioformats
from PIL import Image
import javabridge
javabridge.start_vm(class_path=bioformats.JARS)
from bioformats import ImageReader
import bioformats.omexml as ome
import os
import subprocess

path=('C:\\work\\PyThis\\input')


from os import listdir
from os.path import isfile, join
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

""" notice that even if it just different channels, it may be displayed as different time point or something else
 tempimg2=rawfile.read(c=0,z=0,t=1,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)

"""

xml_string = bioformats.get_omexml_metadata(join(path,onlyfiles[0]))
ome = bioformats.OMEXML(xml_string) # be sure everything is ascii
iome = ome.image(0) # e.g. first image
NumOfC=iome.Pixels.get_SizeC()
NumOfT=iome.Pixels.get_SizeT()
NumOfZ=iome.Pixels.get_SizeZ()
PixelT=iome.Pixels.get_PixelType()

NumOfTReal=int(NumOfT/3);
for x in onlyfiles:
    print (x)
    for t in range(NumOfT):
        print (t)
        timeToFile=int((t-(t%3))/3)
       
        if t%3==0:
            rawfile2 = bioformats.ImageReader('C:\\work\\PyThis\\output\\channel1\\Temp\\Ilastik\\'+x[:-8]+'T%d' %timeToFile +'.tiff')
            tempimg2=rawfile2.read(c=0,z=0,t=0,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
            bioformats.write_image('C:\\work\\PyThis\\output\\channel1\\fused\\'+x[:-8]+'.tiff', tempimg2, PixelT, c=0, z=0, t=timeToFile, size_c=1, size_z=1, size_t=NumOfTReal, channel_names=None)
            rawfile2.close()
            del tempimg2, rawfile2
        else:
            if t%3==1:
                rawfile2 = bioformats.ImageReader('C:\\work\\PyThis\\output\\channel2\\Temp\\Ilastik\\'+x[:-8]+'T%d' %timeToFile +'.tiff')
                tempimg2=rawfile2.read(c=0,z=0,t=0,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
                bioformats.write_image('C:\\work\\PyThis\\output\\channel2\\fused\\'+x[:-8]+'.tiff', tempimg2, PixelT, c=0, z=0, t=timeToFile, size_c=1, size_z=1, size_t=NumOfTReal, channel_names=None)
                rawfile2.close()
            else:
                print (t)
           
 
    
