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
import scipy.ndimage as ndimage
from scipy import stats
import statsmodels.api as sm



myloglevel="ERROR"  # user string argument for logLevel.

javabridge.start_vm(class_path=bioformats.JARS)

rootLoggerName = javabridge.get_static_field("org/slf4j/Logger","ROOT_LOGGER_NAME", "Ljava/lang/String;")
rootLogger = javabridge.static_call("org/slf4j/LoggerFactory","getLogger", "(Ljava/lang/String;)Lorg/slf4j/Logger;", rootLoggerName)
logLevel = javabridge.get_static_field("ch/qos/logback/classic/Level",myloglevel, "Lch/qos/logback/classic/Level;")
javabridge.call(rootLogger, "setLevel", "(Lch/qos/logback/classic/Level;)V", logLevel)



path=('D:\\Zeiss\\wt_vs_S_554_300ng\\S\\input\\')

from os.path import isfile, join
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]


#x=onlyfiles[1]

ResultsS=np.zeros((len(onlyfiles),12)) 
structure_z = np.ones((3, 3, 3), dtype=bool)


for i in range(len(onlyfiles)):
    x=onlyfiles[i]
    xml_string = bioformats.get_omexml_metadata(path+'\\'+x)
    ome = bioformats.OMEXML(xml_string) # be sure everything is ascii
    iome = ome.image(0) # e.g. first image
    NumOfZ=iome.Pixels.get_SizeZ()
    PixelT=iome.Pixels.get_PixelType()
    rawfileiRFP = bioformats.ImageReader(join(path[:-6]+'output\\Channel3_iRFP\\' + x+'f'))
    rawfileNucleusMask = bioformats.ImageReader(join(path[:-6]+'output\\Channel3_iRFP\\Ilastik\\' + x+ 'f'))
    rawfileHubsMask = bioformats.ImageReader(join(path[:-6]+'output\\Cleaned\\Hubs\\' + x+ 'f'))
    rawfileJF = bioformats.ImageReader(join(path[:-6]+'output\\channel1_JF\\' + x+ 'f'))
    rawfileCer = bioformats.ImageReader(join(path[:-6]+'output\\channel2_Cer\\' + x+ 'f'))
    rawfileNocleoliMask = bioformats.ImageReader(join(path[:-6]+'output\\channel2_Cer\\Nucleoli\\' + x+ 'f'))


    print (x)
    MaskNuc=np.zeros((344,344,NumOfZ))
    MaskHubs=np.zeros((344,344,NumOfZ))
    RFPRaw=np.zeros((344,344,NumOfZ))
    JFRaw=np.zeros((344,344,NumOfZ))
    CerRaw=np.zeros((344,344,NumOfZ))
    filled_nucleus=np.zeros((344,344,NumOfZ))
    Nucleoli_mask_nucleoli=np.zeros((344,344,NumOfZ))
    Nucleoli_mask_nucleus=np.zeros((344,344,NumOfZ))
    filled_nucleus_for_nucleoli=np.zeros((344,344,NumOfZ))




    
    for z in range(NumOfZ):  #Ilastik marked the Z slices as T slices, this part loads the images along all z slices
        MaskNuc[:,:,z]=rawfileNucleusMask.read(c=0,z=0,t=z,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        CerRaw[:,:,z]=rawfileCer.read(c=0,z=z,t=0,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        JFRaw[:,:,z]=rawfileJF.read(c=0,z=z,t=0,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        MaskHubs[:,:,z]=rawfileHubsMask.read(c=0,z=0,t=z,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        RFPRaw[:,:,z]=rawfileiRFP.read(c=0,z=z,t=0,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        Nucleoli_mask_nucleus[:,:,z]=rawfileNocleoliMask.read(c=0,z=0,t=z,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)


    rawfileiRFP.close()
    rawfileNucleusMask.close()
    rawfileHubsMask.close()
    rawfileJF.close()
    rawfileCer.close()
    rawfileNocleoliMask.close()
    
    del rawfileiRFP,rawfileNucleusMask,rawfileHubsMask,rawfileJF,rawfileCer,rawfileNocleoliMask
    
    MaskNuc=MaskNuc*(MaskNuc==3)/3 
    Eroded_nuc = ndimage.binary_erosion(MaskNuc, structure=structure_z)
    Eroded_nuc = ndimage.binary_erosion(Eroded_nuc, structure=structure_z)
    Nucleoli_mask_nucleus=Nucleoli_mask_nucleus*(Nucleoli_mask_nucleus==2)/2
    
    for z in range(NumOfZ):
        filled_nucleus[:,:,z] = binary_fill_holes(Eroded_nuc[:,:,z])
        filled_nucleus_for_nucleoli[:,:,z] = binary_fill_holes(Nucleoli_mask_nucleus[:,:,z])

    
    del MaskNuc, Eroded_nuc
    labeled_image, num_features = label(filled_nucleus)
    
    # Find the sizes of each connected component
    component_sizes = np.bincount(labeled_image.ravel())
    # Exclude the zero component (background) and find the largest non-zero component
    non_zero_component_sizes = component_sizes[1:]
    largest_component = non_zero_component_sizes.argmax() + 1
    filled_nucleus=labeled_image*(labeled_image==largest_component)/largest_component
    del labeled_image, component_sizes, num_features,  non_zero_component_sizes, largest_component
    
     
    MaskHubs=MaskHubs*(MaskHubs==2)/2

    Nucleoli_mask_nucleoli=filled_nucleus_for_nucleoli-Nucleoli_mask_nucleus
    Nucleoli_mask_nucleoli=Nucleoli_mask_nucleoli*filled_nucleus
    
    #filter small nucleoli
    labeled_image, num_features = ndimage.label(Nucleoli_mask_nucleoli)
    component_sizes = ndimage.sum(Nucleoli_mask_nucleoli, labeled_image, range(num_features + 1))
    for n in range(len(component_sizes)):
        if component_sizes[n]>350:
            labeled_image[labeled_image==n]=1
        else:
            labeled_image[labeled_image==n]=0
    
    Nucleoli_measure=labeled_image
    del labeled_image, component_sizes ,num_features
    
    Nuc=filled_nucleus-MaskHubs-Nucleoli_mask_nucleoli
    Nuc[Nuc < 0] = 0
    
    JF_in_Nucleoli=JFRaw*Nucleoli_measure
    RFP_in_Nucleoli=RFPRaw*Nucleoli_measure
    
    JF_in_Hubs=JFRaw*MaskHubs
    RFP_in_Hubs=RFPRaw*MaskHubs
    Cer_in_Hubs=CerRaw*MaskHubs
    
    JF_in_Nucleus=JFRaw*Nuc
    RFP_in_Nucleus=RFPRaw*Nuc
    Cer_in_Nucleus=CerRaw*Nuc
    
    TotaliRFPMask=Nuc+MaskHubs+Nucleoli_measure
    TotaliRFPMask=TotaliRFPMask*(TotaliRFPMask==1)/1
    TotaliRFP=TotaliRFPMask*RFPRaw
    TotalJF=TotaliRFPMask*JFRaw
    
    labeled_image, num_of_hubs = label(MaskHubs)
    
    ResultsS[i,0]=num_of_hubs #num of hubs
    ResultsS[i,1]=np.sum(MaskHubs)/num_of_hubs #mean volume
    ResultsS[i,2]=np.sum(JF_in_Hubs)/np.sum(MaskHubs) #mean JF in hubs
    ResultsS[i,3]=np.sum(Cer_in_Hubs)/np.sum(MaskHubs) #mean Cer in hubs
    ResultsS[i,4]=np.sum(RFP_in_Hubs)/np.sum(MaskHubs) #mean RFP in hubs
    
    ResultsS[i,5]=np.sum(JF_in_Nucleus)/np.sum(Nuc) #mean JF in nucleosol
    ResultsS[i,6]=np.sum(Cer_in_Nucleus)/np.sum(Nuc) #mean Cer in nucleosol
    ResultsS[i,7]=np.sum(RFP_in_Nucleus)/np.sum(Nuc) #mean RFP in nucleosol
    
    ResultsS[i,8]=np.sum(JF_in_Nucleoli)/np.sum(Nucleoli_measure) #mean JF in nucleoli
    ResultsS[i,9]=np.sum(RFP_in_Nucleoli)/np.sum(Nucleoli_measure) #mean RFP in nucleolo
        
    ResultsS[i,10]=np.sum(TotalJF)/np.sum(TotaliRFPMask) #total jf
    ResultsS[i,11]=np.sum(TotaliRFP)/np.sum(TotaliRFPMask) #total rfp
    
    
    del filled_nucleus_for_nucleoli, Nucleoli_mask_nucleoli, filled_nucleus, MaskHubs
    del Nuc, JF_in_Nucleoli, RFP_in_Nucleoli, JF_in_Hubs, RFP_in_Hubs, Cer_in_Hubs, Nucleoli_mask_nucleus
    del JF_in_Nucleus, RFP_in_Nucleus, Cer_in_Nucleus, labeled_image, num_of_hubs
    del CerRaw, RFPRaw, JFRaw , TotalJF, TotaliRFP, TotaliRFPMask


  #  plt.imshow(labeled_image[:,:,11])

print('wt')
print('wt')

print('wt')

path=('D:\\Zeiss\\wt_vs_S_554_300ng\\wt\\input\\')

from os.path import isfile, join
onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
Resultswt=np.zeros((len(onlyfiles),12)) 



for i in range(len(onlyfiles)):
    x=onlyfiles[i]
    xml_string = bioformats.get_omexml_metadata(path+'\\'+x)
    ome = bioformats.OMEXML(xml_string) # be sure everything is ascii
    iome = ome.image(0) # e.g. first image
    NumOfZ=iome.Pixels.get_SizeZ()
    PixelT=iome.Pixels.get_PixelType()
    rawfileiRFP = bioformats.ImageReader(join(path[:-6]+'output\\Channel3_iRFP\\' + x+'f'))
    rawfileNucleusMask = bioformats.ImageReader(join(path[:-6]+'output\\Channel3_iRFP\\Ilastik\\' + x+ 'f'))
    rawfileHubsMask = bioformats.ImageReader(join(path[:-6]+'output\\Cleaned\\Hubs\\' + x+ 'f'))
    rawfileJF = bioformats.ImageReader(join(path[:-6]+'output\\channel1_JF\\' + x+ 'f'))
    rawfileCer = bioformats.ImageReader(join(path[:-6]+'output\\channel2_Cer\\' + x+ 'f'))
    rawfileNocleoliMask = bioformats.ImageReader(join(path[:-6]+'output\\channel2_Cer\\Nucleoli\\' + x+ 'f'))


    print (x)
    MaskNuc=np.zeros((344,344,NumOfZ))
    MaskHubs=np.zeros((344,344,NumOfZ))
    RFPRaw=np.zeros((344,344,NumOfZ))
    JFRaw=np.zeros((344,344,NumOfZ))
    CerRaw=np.zeros((344,344,NumOfZ))
    filled_nucleus=np.zeros((344,344,NumOfZ))
    Nucleoli_mask_nucleoli=np.zeros((344,344,NumOfZ))
    Nucleoli_mask_nucleus=np.zeros((344,344,NumOfZ))
    filled_nucleus_for_nucleoli=np.zeros((344,344,NumOfZ))




    
    for z in range(NumOfZ):  #Ilastik marked the Z slices as T slices, this part loads the images along all z slices
        MaskNuc[:,:,z]=rawfileNucleusMask.read(c=0,z=0,t=z,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        CerRaw[:,:,z]=rawfileCer.read(c=0,z=z,t=0,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        JFRaw[:,:,z]=rawfileJF.read(c=0,z=z,t=0,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        MaskHubs[:,:,z]=rawfileHubsMask.read(c=0,z=0,t=z,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        RFPRaw[:,:,z]=rawfileiRFP.read(c=0,z=z,t=0,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)
        Nucleoli_mask_nucleus[:,:,z]=rawfileNocleoliMask.read(c=0,z=0,t=z,series=None,index=None,rescale=False,wants_max_intensity=False,channel_names=None,XYWH=None)


    rawfileiRFP.close()
    rawfileNucleusMask.close()
    rawfileHubsMask.close()
    rawfileJF.close()
    rawfileCer.close()
    rawfileNocleoliMask.close()
    
    del rawfileiRFP,rawfileNucleusMask,rawfileHubsMask,rawfileJF,rawfileCer,rawfileNocleoliMask
    
    MaskNuc=MaskNuc*(MaskNuc==3)/3 
    Eroded_nuc = ndimage.binary_erosion(MaskNuc, structure=structure_z)
    Eroded_nuc = ndimage.binary_erosion(Eroded_nuc, structure=structure_z)
    Nucleoli_mask_nucleus=Nucleoli_mask_nucleus*(Nucleoli_mask_nucleus==2)/2
    
    for z in range(NumOfZ):
        filled_nucleus[:,:,z] = binary_fill_holes(Eroded_nuc[:,:,z])
        filled_nucleus_for_nucleoli[:,:,z] = binary_fill_holes(Nucleoli_mask_nucleus[:,:,z])

    
    del MaskNuc, Eroded_nuc
    labeled_image, num_features = label(filled_nucleus)
    
    # Find the sizes of each connected component
    component_sizes = np.bincount(labeled_image.ravel())
    # Exclude the zero component (background) and find the largest non-zero component
    non_zero_component_sizes = component_sizes[1:]
    largest_component = non_zero_component_sizes.argmax() + 1
    filled_nucleus=labeled_image*(labeled_image==largest_component)/largest_component
    del labeled_image, component_sizes, num_features,  non_zero_component_sizes, largest_component
    
     
    MaskHubs=MaskHubs*(MaskHubs==2)/2

    Nucleoli_mask_nucleoli=filled_nucleus_for_nucleoli-Nucleoli_mask_nucleus
    Nucleoli_mask_nucleoli=Nucleoli_mask_nucleoli*filled_nucleus
    
    #filter small nucleoli
    labeled_image, num_features = ndimage.label(Nucleoli_mask_nucleoli)
    component_sizes = ndimage.sum(Nucleoli_mask_nucleoli, labeled_image, range(num_features + 1))
    for n in range(len(component_sizes)):
        if component_sizes[n]>350:
            labeled_image[labeled_image==n]=1
        else:
            labeled_image[labeled_image==n]=0
    
    Nucleoli_measure=labeled_image
    del labeled_image, component_sizes ,num_features
    
    Nuc=filled_nucleus-MaskHubs-Nucleoli_mask_nucleoli
    Nuc[Nuc < 0] = 0
    
    JF_in_Nucleoli=JFRaw*Nucleoli_measure
    RFP_in_Nucleoli=RFPRaw*Nucleoli_measure
    
    JF_in_Hubs=JFRaw*MaskHubs
    RFP_in_Hubs=RFPRaw*MaskHubs
    Cer_in_Hubs=CerRaw*MaskHubs
    
    JF_in_Nucleus=JFRaw*Nuc
    RFP_in_Nucleus=RFPRaw*Nuc
    Cer_in_Nucleus=CerRaw*Nuc
    TotaliRFPMask=Nuc+MaskHubs+Nucleoli_measure
    TotaliRFPMask=TotaliRFPMask*(TotaliRFPMask==1)/1
    TotaliRFP=TotaliRFPMask*RFPRaw
    TotalJF=TotaliRFPMask*JFRaw


    labeled_image, num_of_hubs = label(MaskHubs)
    
    Resultswt[i,0]=num_of_hubs #num of hubs
    Resultswt[i,1]=np.sum(MaskHubs)/num_of_hubs #mean volume
    Resultswt[i,2]=np.sum(JF_in_Hubs)/np.sum(MaskHubs) #mean JF in hubs
    Resultswt[i,3]=np.sum(Cer_in_Hubs)/np.sum(MaskHubs) #mean Cer in hubs
    Resultswt[i,4]=np.sum(RFP_in_Hubs)/np.sum(MaskHubs) #mean RFP in hubs
    
    Resultswt[i,5]=np.sum(JF_in_Nucleus)/np.sum(Nuc) #mean JF in nucleosol
    Resultswt[i,6]=np.sum(Cer_in_Nucleus)/np.sum(Nuc) #mean Cer in nucleosol
    Resultswt[i,7]=np.sum(RFP_in_Nucleus)/np.sum(Nuc) #mean RFP in nucleosol
    
    Resultswt[i,8]=np.sum(JF_in_Nucleoli)/np.sum(Nucleoli_measure) #mean JF in nucleoli
    Resultswt[i,9]=np.sum(RFP_in_Nucleoli)/np.sum(Nucleoli_measure) #mean RFP in nucleolo
    
    Resultswt[i,10]=np.sum(TotalJF)/np.sum(TotaliRFPMask) 
    Resultswt[i,11]=np.sum(TotaliRFP)/np.sum(TotaliRFPMask)

    
    del filled_nucleus_for_nucleoli, Nucleoli_mask_nucleoli, filled_nucleus, MaskHubs
    del Nuc, JF_in_Nucleoli, RFP_in_Nucleoli, JF_in_Hubs, RFP_in_Hubs, Cer_in_Hubs, Nucleoli_mask_nucleus
    del JF_in_Nucleus, RFP_in_Nucleus, Cer_in_Nucleus, labeled_image, num_of_hubs
    del CerRaw, RFPRaw, JFRaw , TotalJF, TotaliRFP, TotaliRFPMask


X_wt_RFP_nuc=Resultswt[:,7]
X_S_RFP_nuc=ResultsS[:,7]

Y_wt_Vol=Resultswt[:,1]*0.003284566056
Y_S_Vol=ResultsS[:,1]*0.003284566056

X_wt_RFP_nuc_capped=X_wt_RFP_nuc[Y_wt_Vol<10]
Y_wt_Vol_capped=Y_wt_Vol[Y_wt_Vol<10]

X_wt_RFP_nuc_sm= sm.add_constant(X_wt_RFP_nuc)
X_S_RFP_nuc_sm= sm.add_constant(X_S_RFP_nuc)
X_wt_RFP_nuc_capped_sm= sm.add_constant(X_wt_RFP_nuc_capped)

model_wt = sm.OLS(Y_wt_Vol, X_wt_RFP_nuc_sm).fit()
model_S = sm.OLS(Y_S_Vol, X_S_RFP_nuc_sm).fit()
model_wt_capped = sm.OLS(Y_wt_Vol_capped, X_wt_RFP_nuc_capped_sm).fit()


interceptwt, slopewt = model_wt.params
interceptS, slopeS = model_S.params
interceptCapped, slopeCapped = model_wt_capped.params


wt_SEslope = model_wt.bse[1]
S_SEslope = model_S.bse[1]
Capped_SEslope = model_wt_capped.bse[1]


# Calculate the difference in slopes and standard error of the difference
slope_diff = slopewt - slopeS
se_diff = np.sqrt(wt_SEslope**2 + S_SEslope**2)

# Calculate the t-statistic
t_stat = slope_diff / se_diff

# Degrees of freedom
df = len(Y_wt_Vol) + len(Y_S_Vol) - 4  # n1 + n2 - 2*number_of_params

# Calculate the p-value
p_value_wtVsS = 2 * (1 - stats.t.cdf(np.abs(t_stat), df))



# Calculate the difference in slopes and standard error of the difference
slope_diff = slopeCapped - slopeS
se_diff = np.sqrt(Capped_SEslope**2 + S_SEslope**2)

# Calculate the t-statistic
t_stat = slope_diff / se_diff

# Degrees of freedom
df = len(Y_wt_Vol_capped) + len(Y_S_Vol) - 4  # n1 + n2 - 2*number_of_params

# Calculate the p-value
p_value_CappedVsS = 2 * (1 - stats.t.cdf(np.abs(t_stat), df))

#compare the intercept

model_S = sm.OLS(Y_S_Vol, X_S_RFP_nuc_sm).fit()
model_wt_capped = sm.OLS(Y_wt_Vol_capped, X_wt_RFP_nuc_capped_sm).fit()

se_intercept_S = model_S.bse[0]
se_intercept_capped = model_wt_capped.bse[0]

interceptS, slopeS = model_S.params
interceptCapped, slopeCapped = model_wt_capped.params


intercept_diff = interceptCapped - interceptS
se_diff = np.sqrt(se_intercept_capped**2 + se_intercept_S**2)

# Calculate the t-statistic
t_stat = intercept_diff / se_diff

# Degrees of freedom
df = len(Y_wt_Vol_capped) + len(Y_S_Vol) - 4  # n1 + n2 - 2*number_of_params

# Calculate the p-value using the t-distribution
p_value_intercept_capped_vs_S = 2 * (1 - stats.t.cdf(np.abs(t_stat), df))
p_value_intercept_capped_vs_S





#compare the intercept_shifted

X_S_RFP_nuc_sm_50= sm.add_constant(X_S_RFP_nuc-50)
X_wt_RFP_nuc_capped_sm_50= sm.add_constant(X_wt_RFP_nuc_capped-50)

model_S_50 = sm.OLS(Y_S_Vol, X_S_RFP_nuc_sm_50).fit()
model_wt_capped_50 = sm.OLS(Y_wt_Vol_capped, X_wt_RFP_nuc_capped_sm_50).fit()

se_intercept_S = model_S_50.bse[0]
se_intercept_capped = model_wt_capped_50.bse[0]

interceptS_50, slopeS_50 = model_S_50.params
interceptCapped_50, slopeCapped_50 = model_wt_capped_50.params


intercept_diff = interceptCapped_50 - interceptS_50
se_diff = np.sqrt(se_intercept_capped**2 + se_intercept_S**2)

interceptCapped_50 - interceptS_50



# Calculate the t-statistic
t_stat = intercept_diff / se_diff

# Degrees of freedom
df = len(Y_wt_Vol_capped) + len(Y_S_Vol) - 4  # n1 + n2 - 2*number_of_params

# Calculate the p-value using the t-distribution
p_value_intercept_capped_vs_S = 2 * (1 - stats.t.cdf(np.abs(t_stat), df))
p_value_intercept_capped_vs_S




# Get the summary of the model
#print(model_S_20.summary())

# Extract the confidence intervals for the coefficients
#conf_int = model_wt.conf_int(alpha=0.05)  # 95% confidence intervals
#print(conf_int)