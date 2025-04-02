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
projectLoc=('C:\\work\\PyThis\\MyProject.ilp')

from os import listdir
from os.path import isfile, join
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
ilastikfolder=('C:\\Program Files\\ilastik-1.4.0b20-gpu\\')

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
    rawfile = bioformats.ImageReader(join(path,x))
    print (x)
    for t in range(NumOfT):
        if t%6==0:
            print(t)
        
        tempimg=rawfile.read(c=0,z=0,t=t,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        timeToFile=int((t-(t%3))/3)
       
        if t%3==1:
            bioformats.write_image('C:\\work\\PyThis\\output\\channel2\\'+x[:-8]+'.tiff', tempimg, PixelT, c=0, z=0, t=timeToFile, size_c=1, size_z=1, size_t=NumOfTReal, channel_names=None)
        #    bioformats.write_image('C:\\work\\PyThis\\output\\channel2\\Temp\\'+x[:-8]+'T%d' %timeToFile +'.tiff', tempimg, PixelT, c=0, z=0, t=0, size_c=1, size_z=1, size_t=1, channel_names=None)
        else:
            if t%3==2:
                bioformats.write_image('C:\\work\\PyThis\\output\\channel3\\'+x[:-8]+'.tiff', tempimg, PixelT, c=0, z=0, t=timeToFile, size_c=1, size_z=1, size_t=NumOfTReal, channel_names=None)
               # bioformats.write_image('C:\\work\\PyThis\\output\\channel3\\Temp\\'+x[:-8]+'T%d' %timeToFile +'.tiff', tempimg, PixelT, c=0, z=0, t=0, size_c=1, size_z=1, size_t=1, channel_names=None)
            else:
                bioformats.write_image('C:\\work\\PyThis\\output\\channel1\\'+x[:-8]+'.tiff', tempimg, PixelT, c=0, z=0, t=timeToFile, size_c=1, size_z=1, size_t=NumOfTReal, channel_names=None)
            #    bioformats.write_image('C:\\work\\PyThis\\output\\channel1\\Temp\\'+x[:-8]+'T%d' %timeToFile +'.tiff', tempimg, PixelT, c=0, z=0, t=0, size_c=1, size_z=1, size_t=1, channel_names=None)
          
        del tempimg, 
    rawfile.close()
    del rawfile
  
    
#subprocess.run('"'+ ilastikfolder +'run-ilastik.bat" --headless --export_source="Simple Segmentation" --project=' + projectLoc + '--output_filename_format=/IlasPy/{nickname}.tiff' +  'C:\\work\\PyThis\\output\\channel1\\'+x[12:-8]+'.tiff')         



"""
xml_string = bioformats.get_omexml_metadata(join(path,x))
ome = bioformats.OMEXML(xml_string) # be sure everything is ascii
print(ome.image_count)

iome = ome.image(0) # e.g. first image
print (iome.get_Name())

    # get pixel meta data
print (iome.Pixels.get_DimensionOrder())
print (iome.Pixels.get_PixelType())
print (iome.Pixels.get_SizeX())
print (iome.Pixels.get_SizeY())
print (iome.Pixels.get_SizeZ())
print (iome.Pixels.get_SizeT())
print (iome.Pixels.get_SizeC())
print (iome.Pixels.DimensionOrder)
"""
