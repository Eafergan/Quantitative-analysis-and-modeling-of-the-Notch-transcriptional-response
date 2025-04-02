%reset -f
import numpy as np
import matplotlib.pyplot as plt

RBPJmns = np.load('Results_RBPJmns.npy')
RBPJpls = np.load('Results_RBPJpls.npy')
DAPT = np.load('Results_DAPT.npy')

Scatter_RBPJmns = np.load('Scatter_RBPJmns.npy')
Scatter_RBPJpls = np.load('Scatter_RBPJpls.npy')
Scatter_DAPT = np.load('Scatter_DAPT.npy')


Scatter_RBPJmns = Scatter_RBPJmns[0:1499,:]
Scatter_RBPJpls = Scatter_RBPJpls[0:1039,:]
Scatter_DAPT = Scatter_DAPT[0:1241,:]

Norm_Scatter_RBPJmns = np.load('Norm_Scatter_RBPJmns.npy')
Norm_Scatter_RBPJpls = np.load('Norm_Scatter_RBPJpls.npy')
Norm_Scatter_DAPT = np.load('Norm_Scatter_DAPT.npy')


Norm_Scatter_RBPJmns = Norm_Scatter_RBPJmns[0:1499,:]
Norm_Scatter_RBPJpls = Norm_Scatter_RBPJpls[0:1039,:]
Norm_Scatter_DAPT = Norm_Scatter_DAPT[0:1241,:]



# in the result ariable 0:size 1:JFsignal 2:CerSignal 3: num of hubs 4: average nuclear JF sig 5: average nuclear Cer Sig
#Results[:,0,0:2] are averages Results[:,-1,0:2] are st dev
# 0 is DAPT 1 is RBPJ- 2 is RBPJ+

NumOfHubsMat=np.zeros((120,3))
VolOfHubs=np.zeros((120,3))
Norm_Dev_VolOfHubs=np.zeros((120,3))

precentageInHubsJF=np.zeros((120,3))
precentageInHubsCer=np.zeros((120,3))


JF_sig=np.zeros((120,3))
Norm_JF_sig=np.zeros((120,3))
Norm_Dev_JF=np.zeros((120,3))

Cer_sig=np.zeros((120,3))
Norm_Cer_sig=np.zeros((120,3))
Norm_Dev_cer=np.zeros((120,3))

nuclearJF=np.zeros((120,3))
nuclearCer=np.zeros((120,3))

nucleoliJF=np.zeros((120,3))
nucleoliCer=np.zeros((120,3))

eachHub_JF_over_Cer=np.zeros((120,3))


for i in range(len(DAPT)):
    NumOfHub=int(np.max(DAPT[i,:,3]))
    NumOfHubsMat[i,0]=NumOfHub
    nuclearJF[i,0]=DAPT[i,0,4]
    nuclearCer[i,0]=DAPT[i,0,5]
    
    nucleoliJF[i,0]=DAPT[i,1,4]
    nucleoliCer[i,0]=DAPT[i,1,5]
    
    precentageInHubsJF[i,0]=DAPT[i,2,4]
    precentageInHubsCer[i,0]=DAPT[i,2,5]

    eachHub_JF_over_Cer[i,0]=np.mean (   (DAPT[i,1:NumOfHub+1,1]/(DAPT[i,0,4]))/(DAPT[i,1:NumOfHub+1,2]/(DAPT[i,0,5]))  )
    
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
    
    nucleoliJF[i,1]=RBPJmns[i,1,4]
    nucleoliCer[i,1]=RBPJmns[i,1,5]
    
    precentageInHubsJF[i,1]=RBPJmns[i,2,4]
    precentageInHubsCer[i,1]=RBPJmns[i,2,5]
    
    eachHub_JF_over_Cer[i,1]=np.mean (   (RBPJmns[i,1:NumOfHub+1,1]/(RBPJmns[i,0,4]))/(RBPJmns[i,1:NumOfHub+1,2]/(RBPJmns[i,0,5]))  )

    
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
    
    nucleoliJF[i,2]=RBPJpls[i,1,4]
    nucleoliCer[i,2]=RBPJpls[i,1,5]
    
    precentageInHubsJF[i,2]=RBPJpls[i,2,4]
    precentageInHubsCer[i,2]=RBPJpls[i,2,5]
    
    eachHub_JF_over_Cer[i,2]=np.mean (   (RBPJpls[i,1:NumOfHub+1,1]/(RBPJpls[i,0,4]))/(RBPJpls[i,1:NumOfHub+1,2]/(RBPJpls[i,0,5]))  )
    
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


colorvec=['red','blue','magenta','green']
samplevec=['Dll1-FC+, DAPT+','Dll1-FC+','Dll1-FC-','Dll1-FC+, Senexin+']

volOfVoxelInMicron=0.003284566056

#X is signal
plt.rcParams.update({'font.size': 10})
plt.scatter(Scatter_RBPJmns[:,1],Scatter_RBPJmns[:,0]*volOfVoxelInMicron,s=2, linestyle='dashed', color=colorvec[2])
#plt.ylim((0,4000))
#plt.xlim((0,15))
#plt.legend(loc=1,prop={'size': 12})
plt.show()

plt.rcParams.update({'font.size': 10})
plt.scatter(Scatter_RBPJmns[:,0]*volOfVoxelInMicron, Scatter_RBPJmns[:,1],s=2, linestyle='dashed', color=colorvec[2])
#plt.ylim((0,4000))
#plt.xlim((0,15))
#plt.legend(loc=1,prop={'size': 12})
plt.show()

plt.rcParams.update({'font.size': 10})
plt.scatter(Scatter_RBPJmns[:,0]*volOfVoxelInMicron, Norm_Scatter_RBPJmns[:,1],s=2, linestyle='dashed', color=colorvec[2])
#plt.ylim((0,4000))
#plt.xlim((0,15))
#plt.legend(loc=1,prop={'size': 12})
plt.show()


#X is Vol
plt.rcParams.update({'font.size': 10})
plt.scatter(Norm_Scatter_RBPJmns[:,0]*volOfVoxelInMicron,Norm_Scatter_RBPJmns[:,2],s=2, linestyle='dashed', color=colorvec[2])
#plt.ylim((0,4000))
#plt.xlim((0,15))
#plt.legend(loc=1,prop={'size': 12})
plt.show()

plt.figure(figsize=(10, 6))
plt.boxplot(Scatter_RBPJmns[:,0]*volOfVoxelInMicron, vert=True, patch_artist=True, boxprops=dict(facecolor='lightblue'))
plt.title('Box Plot of Sizes')
plt.xlabel('Size')
plt.ylabel('Distribution')
plt.grid(True)
plt.show()

VolumeVec=Scatter_RBPJmns[:,0]*volOfVoxelInMicron
Volume_meanOfEachCell=VolOfHubs[:,1]*volOfVoxelInMicron
Norm_dev_of_vol=Norm_Dev_VolOfHubs[:,1]

np.mean(VolumeVec)

plt.hist(Norm_dev_of_vol, bins=30,density=True, color='lightblue', edgecolor='black')
plt.show()


Q1 = np.percentile(VolumeVec, 25)
Q3 = np.percentile(VolumeVec, 75)

# Calculate IQR
IQR = Q3 - Q1
upper_threshold = Q3 + 3 * IQR

plt.hist(VolumeVec[VolumeVec<upper_threshold], bins=30,density=True, color='lightblue', edgecolor='black')
plt.show()

plt.boxplot(VolumeVec[VolumeVec<upper_threshold], vert=True, patch_artist=True, boxprops=dict(facecolor='lightblue'))
plt.show()

JFVec=Norm_Scatter_RBPJmns[:,1]

Q1 = np.percentile(VolumeVec, 25)
Q3 = np.percentile(VolumeVec, 75)

# Calculate IQR
IQR = Q3 - Q1
upper_threshold = Q3 + 1.5 * IQR

plt.hist(JFVec[JFVec<upper_threshold], bins=30,density=True, color='lightblue', edgecolor='black')
plt.show()

CerVec=Norm_Scatter_RBPJmns[:,2]

Q1 = np.percentile(VolumeVec, 25)
Q3 = np.percentile(VolumeVec, 75)

# Calculate IQR
IQR = Q3 - Q1
upper_threshold = Q3 + 1.5 * IQR

plt.hist(CerVec[CerVec<upper_threshold], bins=30,density=True, color='lightblue', edgecolor='black')
plt.show()