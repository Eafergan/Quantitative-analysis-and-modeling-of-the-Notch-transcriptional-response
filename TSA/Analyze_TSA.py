%reset -f
import numpy as np
import matplotlib.pyplot as plt

TSA = np.load('Results_TSA.npy')


Scatter_TSA = np.load('Scatter_TSA.npy')



Scatter_TSA = Scatter_TSA[0:1499,:]


Norm_Scatter_TSA = np.load('Norm_Scatter_TSA.npy')



Norm_Scatter_TSA = Norm_Scatter_TSA[0:1499,:]



# in the result ariable 0:size 1:JFsignal 2:CerSignal 3: num of hubs 4: average nuclear JF sig 5: average nuclear Cer Sig
#Results[:,0,0:2] are averages Results[:,-1,0:2] are st dev
# 0 is DAPT 1 is RBPJ- 2 is RBPJ+

NumOfHubsMat=np.zeros((120))
VolOfHubs=np.zeros((120))
Norm_Dev_VolOfHubs=np.zeros((120))

precentageInHubsJF=np.zeros((120))
precentageInHubsCer=np.zeros((120))


VolOfCell=np.zeros((120))
VolOfCell_wo_nuc=np.zeros((120))

JF_sig=np.zeros((120))
Norm_JF_sig=np.zeros((120))
Norm_Dev_JF=np.zeros((120))

Cer_sig=np.zeros((120))
Norm_Cer_sig=np.zeros((120))
Norm_Dev_cer=np.zeros((120))

nuclearJF=np.zeros((120))
nuclearCer=np.zeros((120))

nucleoliJF=np.zeros((120))
nucleoliCer=np.zeros((120))

eachHub_JF_over_Cer=np.zeros((120))


for i in range(len(TSA)):
    NumOfHub=int(np.max(TSA[i,:,3]))
    NumOfHubsMat[i]=NumOfHub
    nuclearJF[i]=TSA[i,0,4]
    nuclearCer[i]=TSA[i,0,5]
    VolOfCell[i]=TSA[i,3,4]#vol of cell for control
    VolOfCell_wo_nuc[i]=TSA[i,3,5]#vol of cellw woo nocleo for control
    nucleoliJF[i]=TSA[i,1,4]
    nucleoliCer[i]=TSA[i,1,5]
    
    precentageInHubsJF[i]=TSA[i,2,4]
    precentageInHubsCer[i]=TSA[i,2,5]

    eachHub_JF_over_Cer[i]=np.mean (   (TSA[i,1:NumOfHub+1,1]/(TSA[i,0,4]))/(TSA[i,1:NumOfHub+1,2]/(TSA[i,0,5]))  )
    
    VolOfHubs[i]=np.mean(TSA[i,1:NumOfHub+1,0])
    Norm_Dev_VolOfHubs[i]=np.std(TSA[i,1:NumOfHub+1,0])/VolOfHubs[i]
    
    JF_sig[i]=np.mean(TSA[i,1:NumOfHub+1,1])
    Norm_JF_sig[i]=JF_sig[i]/nuclearJF[i]
    Norm_Dev_JF[i]=np.std(TSA[i,1:NumOfHub+1,1])/JF_sig[i] #norm to hub signal
    
    Cer_sig[i]=np.mean(TSA[i,1:NumOfHub+1,2])
    Norm_Cer_sig[i]=Cer_sig[i]/nuclearCer[i]
    Norm_Dev_cer[i]=np.std(TSA[i,1:NumOfHub+1,2])/Cer_sig[i]
    
    

TSA_JF_over_cer=np.zeros((2000,2))
for i in range(2000): #this bootstrap calculate how much weaker the signal after TSA
#JFtsa-:mean= 1.870979828 sterr=0.015846335; JFtsa+:mean=1.792953182, sterr=0.025376414
#CerTSA-: mean= 1.557155345 sterr= 0.016148793; CerTSA+ mean=1.377661364 sterr=0.015695493
    randJFmns=np.random.normal(1.870979828,0.015846335)
    randJFpls=np.random.normal(1.792953182,0.025376414)
    randCermns=np.random.normal(1.557155345,0.016148793)
    randCerpls=np.random.normal(1.377661364,0.015695493)
    TSA_JF_over_cer[i,0]=(randJFmns-randJFpls)/randJFmns
    TSA_JF_over_cer[i,1] =(randCermns-randCerpls)/randCermns

meanJF_how_weak=np.mean(TSA_JF_over_cer[:,0])
meanCer_how_weak=np.mean(TSA_JF_over_cer[:,1])

minJF_how_weak=np.percentile(TSA_JF_over_cer[:,0],2.5)
maxJF_how_weak=np.percentile(TSA_JF_over_cer[:,0],97.5)

minCer_how_weak=np.percentile(TSA_JF_over_cer[:,1],2.5)
maxCer_how_weak=np.percentile(TSA_JF_over_cer[:,1],97.5)

PlusMinusSignJF=(maxJF_how_weak-minJF_how_weak)/2
PlusMinusSignCer=(maxCer_how_weak-minCer_how_weak)/2
