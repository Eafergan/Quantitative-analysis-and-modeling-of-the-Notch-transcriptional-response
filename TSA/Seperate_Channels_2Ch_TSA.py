%reset -f

import bioformats
import javabridge
javabridge.start_vm(class_path=bioformats.JARS)
import os
from os import listdir
from os.path import isfile, join

myloglevel="ERROR"  # user string argument for logLevel.

javabridge.start_vm(class_path=bioformats.JARS)

rootLoggerName = javabridge.get_static_field("org/slf4j/Logger","ROOT_LOGGER_NAME", "Ljava/lang/String;")
rootLogger = javabridge.static_call("org/slf4j/LoggerFactory","getLogger", "(Ljava/lang/String;)Lorg/slf4j/Logger;", rootLoggerName)
logLevel = javabridge.get_static_field("ch/qos/logback/classic/Level",myloglevel, "Lch/qos/logback/classic/Level;")
javabridge.call(rootLogger, "setLevel", "(Lch/qos/logback/classic/Level;)V", logLevel)



path=('D:\\Zeiss\\TSA2uM\\input\\')

from os.path import isfile, join
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]


for x in onlyfiles:
    xml_string = bioformats.get_omexml_metadata(path+'\\'+x)
    ome = bioformats.OMEXML(xml_string) # be sure everything is ascii
    iome = ome.image(0) # e.g. first image
    NumOfZ=iome.Pixels.get_SizeZ()
    PixelT=iome.Pixels.get_PixelType()
    
    rawfile = bioformats.ImageReader(join(path,x))
    print (x)
    for z in range(NumOfZ):
        tempimg1=rawfile.read(c=0,z=z,t=0,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        bioformats.write_image(path[0:-6] + 'output\\channel1_JF\\'+x[:-4]+'.tiff', tempimg1, PixelT, c=0, z=z, t=0, size_c=1, size_z=NumOfZ, size_t=1, channel_names=None)
        tempimg2=rawfile.read(c=1,z=z,t=0,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        bioformats.write_image(path[0:-6] + 'output\\channel2_Cer\\'+x[:-4]+'.tiff', tempimg2, PixelT, c=0, z=z, t=0, size_c=1, size_z=NumOfZ, size_t=1, channel_names=None)          
        del tempimg1, tempimg2
    rawfile.close()
    del rawfile , xml_string, ome, iome, NumOfZ, PixelT
del x


