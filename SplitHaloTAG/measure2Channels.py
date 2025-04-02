%reset -f

from matplotlib import pyplot as plt
import numpy as np
import bioformats
from PIL import Image
import cv2
from skimage.morphology import (erosion, dilation, opening, closing, disk)
import javabridge
javabridge.start_vm(class_path=bioformats.JARS)
from bioformats import ImageReader
import bioformats.omexml as ome
from os import listdir
import xlsxwriter
from datetime import (date ,datetime)

now = datetime.now()
dt_string = now.strftime("%d_%m_%Y_%H_%M_%S")


workbook = xlsxwriter.Workbook(dt_string+'.xlsx')
worksheet = workbook.add_worksheet()
path=('C:\\work\\PyThis\\output\\channel1')
pathread=('C:\\work\\PyThis\\output\\')

from os.path import isfile, join
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]

nameoffile=onlyfiles[0] # just inorder to have first estimaion of times

xml_string = bioformats.get_omexml_metadata(path+'\\'+nameoffile)
ome = bioformats.OMEXML(xml_string) # be sure everything is ascii
iome = ome.image(0) # e.g. first image
NumOfT=iome.Pixels.get_SizeT()
PixelT=iome.Pixels.get_PixelType()

#NumOfT=13

del nameoffile

nameoffile=  onlyfiles[0]

for x in range(48):
    worksheet.write(8*2  + 2, x+1, x) 

for x in range(1): 
    worksheet.write(0+6*x, 0, nameoffile[:-5])  
    worksheet.write(1+6*x, 0, 'Background') 
    worksheet.write(2+6*x,0 , 'Halo')
    worksheet.write(3+6*x,0 , 'Signal To Background')
    worksheet.write(4+6*x,0 , 'Normalized Signal')
    worksheet.write(7*2 + x + 1, 0, nameoffile[:-5]) 
    
    worksheet.write(8*2 + x + 3, 0, nameoffile[:-5]) 
    print(nameoffile)
    
    
    
    
    
    for t in range(NumOfT):
        rawfile = bioformats.ImageReader(path+'\\fused\\'+nameoffile)
        tempmask=rawfile.read(c=0,z=0,t=t,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        
        rawfile1 = bioformats.ImageReader(pathread+'channel1\\'+nameoffile)
        tempmeasure1=rawfile1.read(c=0,z=0,t=t,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)

        
        rawfile2 = bioformats.ImageReader(pathread+'channel2\\'+nameoffile)
        tempmeasure2=rawfile2.read(c=0,z=0,t=t,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        
        mask=tempmask*(tempmask==2)/2
        bgmask=tempmask*(tempmask==1)
        measure1=tempmeasure1*mask
        measure2=tempmeasure2*mask
        bg1=bgmask*tempmeasure1
        bg2=bgmask*tempmeasure2
        
        averageSig1=np.sum(np.float64(measure1))/np.sum(np.float64(mask))
        averageBg1=np.sum(np.float64(bg1))/np.sum(np.float64(bgmask))
        
        averageSig2=np.sum(np.float64(measure2))/np.sum(np.float64(mask))
        averageBg2=np.sum(np.float64(bg2))/np.sum(np.float64(bgmask))
        
      
        
        worksheet.write(2+7*0, t+1, averageSig1)  #this write the value to the excel
        worksheet.write(1+7*0, t+1, averageBg1)          #this write the time
        worksheet.write(3+7*0, t+1, averageSig1-averageBg1) 
        worksheet.write(4+7*0, t+1, (averageSig1-averageBg1)/averageBg1)
        worksheet.write(7*2 + 0 + 1, t+1, averageSig1-averageBg1) 
        worksheet.write(8*2 + 0 + 3, t+1, (averageSig1-averageBg1)/averageBg1) 

        worksheet.write(2+7*1, t+1, averageSig2)  #this write the value to the excel
        worksheet.write(1+7*1, t+1, averageBg2)          #this write the time
        worksheet.write(3+7*1, t+1, averageSig2-averageBg2) 
        worksheet.write(4+7*1, t+1, (averageSig2-averageBg2)/averageBg2)
        worksheet.write(7*2 + 1 + 1, t+1, averageSig2-averageBg2) 
        worksheet.write(8*2 + 1 + 3, t+1, (averageSig2-averageBg2)/averageBg2) 

        print(t)
        del rawfile,rawfile1, rawfile2, tempmask, tempmeasure1, tempmeasure2, mask, bgmask, measure1, measure2, bg1, bg2, averageSig1, averageSig2,  averageBg1, averageBg2

   


workbook.close()
"""
xml_string = bioformats.get_omexml_metadata(path+'\\ilast\\'+nameoffile)
ome = bioformats.OMEXML(xml_string) # be sure everything is ascii
iome = ome.image(0) # e.g. first image
NumOfT=iome.Pixels.get_SizeT()
NumOfZ=iome.Pixels.get_SizeZ()
NumOfC=iome.Pixels.get_SizeC()
PixelT=iome.Pixels.get_PixelType()
"""
""" use this to see image
img=Image.fromarray(opened*100)
img.show()
"""

#bioformats.write_image('C:\\work\\pythis\\output\\a3.tif', bg,  'uint16', c=0, z=0, t=0, size_c=1, size_z=1, size_t=1, channel_names=None)
 