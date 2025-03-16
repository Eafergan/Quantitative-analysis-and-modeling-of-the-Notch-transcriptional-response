%reset -f
import numpy as np
import matplotlib.pyplot as plt

RBPJmns = np.load('Results_RBPJmns.npy')
RBPJpls = np.load('Results_RBPJpls.npy')
DAPT = np.load('Results_DAPT.npy')

# in the result ariable 0:size 1:JFsignal 2:CerSignal 3: num of hubs 4: average nuclear JF sig 5: average nuclear Cer Sig
#Results[:,0,0:2] are averages Results[:,-1,0:2] are st dev
# 0 is DAPT 1 is RBPJ- 2 is RBPJ+

NumOfHubsMat=np.zeros((120,3))
VolOfHubs=np.zeros((120,3))
Norm_Dev_VolOfHubs=np.zeros((120,3))

JF_sig=np.zeros((120,3))
Norm_JF_sig=np.zeros((120,3))
Norm_Dev_JF=np.zeros((120,3))

Cer_sig=np.zeros((120,3))
Norm_Cer_sig=np.zeros((120,3))
Norm_Dev_cer=np.zeros((120,3))

nuclearJF=np.zeros((120,3))
nuclearCer=np.zeros((120,3))

for i in range(len(DAPT)):
    NumOfHub=int(np.max(DAPT[i,:,3]))
    NumOfHubsMat[i,0]=NumOfHub
    nuclearJF[i,0]=DAPT[i,0,4]
    nuclearCer[i,0]=DAPT[i,0,5]
    
    VolOfHubs[i,0]=np.mean(DAPT[i,1:NumOfHub+1,0])
    Norm_Dev_VolOfHubs[i,0]=np.std(DAPT[i,1:NumOfHub+1,0])/VolOfHubs[i,0]
    
    JF_sig[i,0]=np.mean(DAPT[i,1:NumOfHub+1,1])
    Norm_JF_sig[i,0]=JF_sig[i,0]/nuclearJF[i,0]
    Norm_Dev_JF[i,0]=np.std(DAPT[i,1:NumOfHub+1,1])/JF_sig[i,0] #norm to hub signal
    
    Cer_sig[i,0]=np.mean(DAPT[i,1:NumOfHub+1,2])
    Norm_Cer_sig[i,0]=Cer_sig[i,0]/nuclearCer[i,0]
    Norm_Dev_cer[i,0]=np.std(DAPT[i,1:NumOfHub+1,2])/Cer_sig[i,0]
    

for i in range(len(RBPJmns)):
    NumOfHub=int(np.max(RBPJmns[i,:,3]))
    NumOfHubsMat[i,1]=NumOfHub
    nuclearJF[i,1]=RBPJmns[i,0,4]
    nuclearCer[i,1]=RBPJmns[i,0,5]
    
    VolOfHubs[i,1]=np.mean(RBPJmns[i,1:NumOfHub+1,0])
    Norm_Dev_VolOfHubs[i,1]=np.std(RBPJmns[i,1:NumOfHub+1,0])/VolOfHubs[i,1]
    
    JF_sig[i,1]=np.mean(RBPJmns[i,1:NumOfHub+1,1])
    Norm_JF_sig[i,1]=JF_sig[i,1]/nuclearJF[i,1]
    Norm_Dev_JF[i,1]=np.std(RBPJmns[i,1:NumOfHub+1,1])/JF_sig[i,1] #norm to hub signal
    
    Cer_sig[i,1]=np.mean(RBPJmns[i,1:NumOfHub+1,2])
    Norm_Cer_sig[i,1]=Cer_sig[i,1]/nuclearCer[i,1]
    Norm_Dev_cer[i,1]=np.std(RBPJmns[i,1:NumOfHub+1,2])/Cer_sig[i,1]

    
for i in range(len(RBPJpls)):
    NumOfHub=int(np.max(RBPJpls[i,:,3]))
    NumOfHubsMat[i,2]=NumOfHub
    nuclearJF[i,2]=RBPJpls[i,0,4]
    nuclearCer[i,2]=RBPJpls[i,0,5]
    
    VolOfHubs[i,2]=np.mean(RBPJpls[i,1:NumOfHub+1,0])
    Norm_Dev_VolOfHubs[i,2]=np.std(RBPJpls[i,1:NumOfHub+1,0])/VolOfHubs[i,2]
    
    JF_sig[i,2]=np.mean(RBPJpls[i,1:NumOfHub+1,1])
    Norm_JF_sig[i,2]=JF_sig[i,2]/nuclearJF[i,2]
    Norm_Dev_JF[i,2]=np.std(RBPJpls[i,1:NumOfHub+1,1])/JF_sig[i,2] #norm to hub signal
    
    Cer_sig[i,2]=np.mean(RBPJpls[i,1:NumOfHub+1,2])
    Norm_Cer_sig[i,2]=Cer_sig[i,2]/nuclearCer[i,2]
    Norm_Dev_cer[i,2]=np.std(RBPJpls[i,1:NumOfHub+1,2])/Cer_sig[i,2]    


VolOfHubs=np.nan_to_num(VolOfHubs, nan=0.0)
Norm_Dev_VolOfHubs=np.nan_to_num(Norm_Dev_VolOfHubs, nan=0.0)




