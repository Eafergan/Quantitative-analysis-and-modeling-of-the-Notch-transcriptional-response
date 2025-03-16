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
from os import listdir
from datetime import (date ,datetime)
import h5py
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import pandas as pd


myloglevel="ERROR"  # user string argument for logLevel.

javabridge.start_vm(class_path=bioformats.JARS)

rootLoggerName = javabridge.get_static_field("org/slf4j/Logger","ROOT_LOGGER_NAME", "Ljava/lang/String;")
rootLogger = javabridge.static_call("org/slf4j/LoggerFactory","getLogger", "(Ljava/lang/String;)Lorg/slf4j/Logger;", rootLoggerName)
logLevel = javabridge.get_static_field("ch/qos/logback/classic/Level",myloglevel, "Lch/qos/logback/classic/Level;")
javabridge.call(rootLogger, "setLevel", "(Lch/qos/logback/classic/Level;)V", logLevel)



def exp_decay(x, C1,C2,C3): #C1 is ratio, C2 is Lambda for 1 and C3 is for 2, different than in the PPT
    return  C1*np.exp(-x*C2)+ (1-C1)*np.exp(-x*C3)


def Single_exp_decay(x, C1,C2): 
    return  (1-C2)*np.exp(-x*C1) + C2




path=('D:\\PulseChaseOnDll\\E100nM\\input\\')

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

#the image was cropped in FIJI due to a physical artifact  
  
#the following segment fuse both channels


for x in onlyfiles:
    rawfile1 = bioformats.ImageReader(path[0:-6] + 'output\\channel1\\'+x[:-8]+'.tiff')
    rawfile2 = bioformats.ImageReader(path[0:-6] + 'output\\channel2\\'+x[:-8]+'.tiff')
    print (x)
    for t in range(NumOfTReal):
        print (t)
        tempimg1=rawfile1.read(c=0,z=0,t=t,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        tempimg2=rawfile2.read(c=0,z=0,t=t,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        timeToFile=t
        fused=(tempimg1+tempimg2)*0.5
        bioformats.write_image(path[0:-6] + 'output\\fused\\'+x[:-8]+'.tiff', fused, PixelT, c=0, z=0, t=t, size_c=1, size_z=1, size_t=NumOfTReal, channel_names=None)
        del tempimg1, tempimg2 , fused
    rawfile1.close()
    rawfile2.close()
    del rawfile1,  rawfile2 


resultsCh1=np.zeros((NumOfTReal,len(onlyfiles) )) 
resultsCh2=np.zeros((NumOfTReal,len(onlyfiles) )) 
bgCh1=np.zeros((NumOfTReal,len(onlyfiles) )) 
bgCh2=np.zeros((NumOfTReal,len(onlyfiles) )) 
MaskCh2=np.zeros((NumOfTReal,len(onlyfiles) )) 


for x in range(len(onlyfiles)):
    nameoffile=  onlyfiles[x]  
    print(nameoffile[:-8]) 
    Ilastik = h5py.File(path[0:-6] + 'output\\fused\\Ilastik\\'+nameoffile[:-8]+'.h5')
    Mask = Ilastik['exported_data']
    
    for t in range(NumOfTReal):
        tempmask=Mask[t,:,:,0]
        tempbgmask=Mask[t,:,:,0]+1
        rawfile1 = bioformats.ImageReader( path[0:-6] + 'output\\channel1\\'+nameoffile[:-8]+'.tiff')
        rawfile2 = bioformats.ImageReader( path[0:-6] + 'output\\channel2\\'+nameoffile[:-8]+'.tiff')
        tempmeasure1=rawfile1.read(c=0,z=0,t=t,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        tempmeasure2=rawfile2.read(c=0,z=0,t=t,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)

        tempmask=tempmask*(tempmask==3)/3
        tempbgmask=tempbgmask*(tempbgmask==1)*1
        measureMasked1=tempmeasure1*tempmask
        measureMasked2=tempmeasure2*tempmask

        bgMasked1=tempbgmask*tempmeasure1
        bgMasked2=tempbgmask*tempmeasure2

        Sig90Ch1=measureMasked1[(measureMasked1<np.percentile(measureMasked1[measureMasked1>0], 90)) & (measureMasked1>0)]
        Bg90Ch1=bgMasked1[(bgMasked1<np.percentile(bgMasked1[bgMasked1>0], 90)) & (bgMasked1>0)]

        Sig90Ch2=measureMasked2[(measureMasked2<np.percentile(measureMasked2[measureMasked2>0], 90)) & (measureMasked2>0)]
        Bg90Ch2=bgMasked2[(bgMasked2<np.percentile(bgMasked2[bgMasked2>0], 90)) & (bgMasked2>0)]

        averageSig1=np.sum(np.float64(Sig90Ch1))/np.sum(np.float64(tempmask))
        averageBg1=np.sum(np.float64(Bg90Ch1))/np.sum(np.float64(tempbgmask)) ##maybe change size to np.count_nonzero
        averageSig2=np.sum(np.float64(Sig90Ch2))/np.sum(np.float64(tempmask))
        averageBg2=np.sum(np.float64(Bg90Ch2))/np.sum(np.float64(tempbgmask)) ##maybe change size to np.count_nonzero
 
        resultsCh1[t,x]=averageSig1-averageBg1
        resultsCh2[t,x]=averageSig2-averageBg2
        bgCh1[t,x]=averageBg1
        bgCh2[t,x]=averageBg2
        MaskCh2[t,x]=np.sum(np.float64(tempmask))


        if t==0:
            if x==0:
                bioformats.write_image(path[0:-6] + 'output\\measureCh1.tif', measureMasked1,  'uint16', c=0, z=0, t=0, size_c=1, size_z=1, size_t=1, channel_names=None)
                bioformats.write_image(path[0:-6] + 'output\\bgCh1.tif', bgMasked1,  'uint16', c=0, z=0, t=0, size_c=1, size_z=1, size_t=1, channel_names=None)
                bioformats.write_image(path[0:-6] + 'output\\measureCh2.tif', measureMasked2,  'uint16', c=0, z=0, t=0, size_c=1, size_z=1, size_t=1, channel_names=None)
                bioformats.write_image(path[0:-6] + 'output\\bgCh2.tif', bgMasked2,  'uint16', c=0, z=0, t=0, size_c=1, size_z=1, size_t=1, channel_names=None)
     
  
        if t%6==0:
            print(t)
        del rawfile1, rawfile2, tempmask, tempmeasure1,tempmeasure2,  measureMasked1, measureMasked2, averageSig1, averageBg1, averageSig2, averageBg2, bgMasked1, bgMasked2, Sig90Ch1, Sig90Ch2, Bg90Ch1, Bg90Ch2
    del nameoffile, Mask, Ilastik


javabridge.kill_vm()

NamesVector = [elem[elem.find('__')+2:elem.find('.ome')] for elem in onlyfiles]


### so far the results were "resultsCh1", now the for these results to exponential decay.


df2 = pd.read_excel('D:\\PulseChaseOnDll\\CDE.xlsx',sheet_name='Sheet2')
df2AsArray=df2.to_numpy()
flippedArray=np.swapaxes(df2AsArray,0,1)

samples=['DLL+DAPT','DLL','PBS']
normalized_matrix=np.zeros((NumOfTReal,3 )) 
ErrorMat=np.zeros((NumOfTReal,3 )) 
MaskGrowthMat=np.zeros((NumOfTReal,3 )) 
MaskGrowthErrorMat=np.zeros((NumOfTReal,3 )) 


normalized_matrix[:,:]=flippedArray[1:146,20:23]
ErrorMat[:,:]=flippedArray[1:146,26:29]

MaskGrowthMat[:,:]=flippedArray[1:146,50:53] 
MaskGrowthErrorMat[:,:]=flippedArray[1:146,56:59] 



vec = np.linspace(0, (NumOfTReal-1)*10, num=NumOfTReal)

C1Mat=np.zeros(3)
C2Mat=np.zeros(3)
C3Mat=np.zeros(3)
StErrorC1=np.zeros(3)
StErrorC2=np.zeros(3)
StErrorC3=np.zeros(3)
for x in range(3):
    params, covariance = curve_fit(exp_decay, vec, normalized_matrix[:, x],bounds=([0.6,0,0],[1,1,1]))
    y_pred = exp_decay( vec, *params)
    r_squared = r2_score( normalized_matrix[:, x], y_pred)
    Sterr=np.sqrt(np.diag(covariance))
    C1Mat[x]=params[0]
    C2Mat[x]=params[1]
    C3Mat[x]=params[2]
    StErrorC1[x]=Sterr[0]
    StErrorC2[x]=Sterr[1]
    StErrorC3[x]=Sterr[2]
    print(x)
    del params, covariance, y_pred, r_squared, Sterr
    
params_sin, covariance_sin = curve_fit(Single_exp_decay, vec, normalized_matrix[:, 1],bounds=([0,0],[1,1]) )
C1Mat_sin=params_sin[0]
C2Mat_sin=params_sin[1]


            
for i in [0,1,2]:
    y_pred = exp_decay( vec,  C1Mat[i],C2Mat[i], C3Mat[i])
    r_squared = r2_score( normalized_matrix[ :, i], y_pred)
    plt.plot(vec,y_pred,linestyle='dashed')
    plt.fill_between(vec, normalized_matrix[ :, i]-ErrorMat[ :, i], normalized_matrix[ :, i]+ErrorMat[ :, i],alpha=0.3,label='_nolegend_')
    plt.legend(['Dll+ DAPT+','Dll+','Control'])
plt.ylim((0,1))
plt.show()

#single vs two exponents on Dll sample
y_pred = exp_decay( vec,  C1Mat[1],C2Mat[1], C3Mat[1])
y_pred_single = Single_exp_decay( vec,  C1Mat_sin,C2Mat_sin )
y_long = Single_exp_decay( vec, C2Mat[1],0)
y_Short = Single_exp_decay( vec, C3Mat[1],0)
#['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'] color order
plt.plot(vec,y_pred,linestyle='dashed',color='#ff7f0e',label='Fit of two exponents')
plt.plot(vec,C1Mat[1]* y_long,linestyle='dashed', color='green', label='exponent #1')
plt.plot(vec,(C1Mat[1])+(1-C1Mat[1])*y_Short,linestyle='dashed',label='exponent #2')
plt.plot(vec,y_pred_single,linestyle='dashed', color='red', label='Fit of single exponents')
plt.fill_between(vec, normalized_matrix[ :, 1]-ErrorMat[ :, 1], normalized_matrix[ :, 1]+ErrorMat[ :, 1],alpha=0.3,label='_nolegend_', color='#ff7f0e')
plt.legend()
plt.ylim((0,1))
plt.show()

HalfLifeC2=np.log(2)/C2Mat[:]
HalfLifeC3=np.log(2)/C3Mat[:]


for i in [0,1,2]:
    plt.fill_between(vec, MaskGrowthMat[ :, i]-MaskGrowthErrorMat[ :, i], MaskGrowthMat[ :, i]+MaskGrowthErrorMat[ :, i],alpha=0.3)
    plt.legend(['DLL+DAPT','DLL','PBS'])
plt.show()


## old
first_time_points = resultsCh2[0,:]
normalized_matrix = resultsCh2 / first_time_points[np.newaxis, :]



HalfMat=np.zeros(len(onlyfiles))
R2Mat=np.zeros(len(onlyfiles))
C1Mat=np.zeros(len(onlyfiles))
C2Mat=np.zeros(len(onlyfiles))
C3Mat=np.zeros(len(onlyfiles))


AverageMat=np.zeros((NumOfTReal,len(onlyfiles) )) #0is DAPT 1 is dE 2 is dR 3 is H2B
DeviationMat=np.zeros((NumOfTReal,len(onlyfiles) )) #0is DAPT 1 is dE 2 is dR 3 is H2B


AverageMat[:,0]=np.average(normalized_matrix[:,[0,4,8,12]],axis=1)
AverageMat[:,1]=np.average(normalized_matrix[:,[1,5,9,13]],axis=1)
AverageMat[:,2]=np.average(normalized_matrix[:,[2,6,10,14]],axis=1)

first_time_points_Avg = AverageMat[0,:]
normalized_matrix_Avg = AverageMat / first_time_points_Avg[np.newaxis, :]


DeviationMat[:,0]=np.std(normalized_matrix[:,[0,4,8,12]],axis=1)
DeviationMat[:,1]=np.std(normalized_matrix[:,[1,5,9,13]],axis=1)
DeviationMat[:,2]=np.std(normalized_matrix[:,[2,6,10,14]],axis=1)
DeviationMat[:,3]=np.std(normalized_matrix[:,[3,7,11,15]],axis=1)
ErrorMat=DeviationMat/2




first_time_points_Avg = AverageMat[0,:]
normalized_matrix_Avg = AverageMat / first_time_points_Avg[np.newaxis, :]


R2Mat_Avg=np.zeros(len(onlyfiles))
C1Mat_Avg=np.zeros(len(onlyfiles))
C2Mat_Avg=np.zeros(len(onlyfiles))
C3Mat_Avg=np.zeros(len(onlyfiles))


for x in range(len(onlyfiles)):
    params, covariance = curve_fit(exp_decay, vec, normalized_matrix_Avg[:, x], bounds=([0.6,0,0],[1,1,1]))
    y_pred = exp_decay( vec, *params)
    r_squared = r2_score( normalized_matrix[:, x], y_pred)
    R2Mat_Avg[x]=r_squared
    C1Mat_Avg[x]=params[0]
    C2Mat_Avg[x]=params[1]
    C3Mat_Avg[x]=params[2]
    del params, covariance, y_pred, r_squared
 

#use this for test
fileChosen=2
y_pred = exp_decay( vec,  C1Mat[fileChosen],C2Mat[fileChosen], C3Mat[fileChosen])
r_squared = r2_score( normalized_matrix[ :, fileChosen], y_pred)
plt.plot(range(0,145),normalized_matrix[:,fileChosen])
plt.plot(range(0,145),y_pred)
plt.ylabel('')
plt.show()

#use this for test
fileChosen=1
y_pred = exp_decay( vec,  C1Mat_Avg[fileChosen],C2Mat_Avg[fileChosen], C3Mat_Avg[fileChosen])
r_squared = r2_score( normalized_matrix_Avg[ :, fileChosen], y_pred)
plt.plot(vec,normalized_matrix_Avg[:,fileChosen])
plt.fill_between(vec, normalized_matrix_Avg[ :, fileChosen]-ErrorMat[ :, fileChosen], normalized_matrix_Avg[ :, fileChosen]+ErrorMat[ :, fileChosen])
plt.plot(vec,y_pred)
plt.ylabel('')
plt.show()

for i in [1,3]:
    y_pred = exp_decay( vec,  C1Mat_Avg[i],C2Mat_Avg[i], C3Mat_Avg[i])
    r_squared = r2_score( normalized_matrix_Avg[ :, i], y_pred)
    plt.plot(vec,y_pred,linestyle='dashed')
    plt.fill_between(vec, normalized_matrix_Avg[ :, i]-ErrorMat[ :, i], normalized_matrix_Avg[ :, i]+ErrorMat[ :, i],alpha=0.3,label='_nolegend_')
    plt.legend(['Notch1dE','H2B Control'])
plt.show()


HalfLife=np.log(2)/0.015781


 #C1 is ratio, C2 is Lambda for 1 and C3 is for 2