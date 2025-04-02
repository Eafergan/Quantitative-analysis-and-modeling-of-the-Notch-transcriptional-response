%reset -f


from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.optimize import minimize


def exp_recovary(x, C1,C2,C3): 
    return  C2*(1-np.exp(-x*C1))*np.exp(-C3*x)

def T_half(x, C1, C2, C3, target_value):
    return ((exp_recovary(x, C1, C2, C3)/C2) - target_value)**2

#####

Con_nuc_df = pd.read_excel(r'D:\Zeiss\Frap\FrapXLS.xlsx',sheet_name='Control size of nucFrap')
Con_nuc=Con_nuc_df.to_numpy()

Con_hub_df = pd.read_excel(r'D:\Zeiss\Frap\FrapXLS.xlsx',sheet_name='Control Size Hub')
Con_hub=Con_hub_df.to_numpy()

nuc_df = pd.read_excel(r'D:\Zeiss\Frap\FrapXLS.xlsx',sheet_name='Nucleosol')
nuc=nuc_df.to_numpy()

hub_df = pd.read_excel(r'D:\Zeiss\Frap\FrapXLS.xlsx',sheet_name='hubs')
hub=hub_df.to_numpy()

del Con_nuc_df, Con_hub_df, nuc_df, hub_df



#####

#####
Con_nuc_DAPT=Con_nuc[2:602,0:12]
Con_nuc_RBPJmns=Con_nuc[2:602,14:26]
Con_nuc_RBPJpls=Con_nuc[2:602,28:40]
Con_nuc_NICD=Con_nuc[2:602,42:54]
del Con_nuc

Con_hub_DAPT=Con_hub[2:602,0:12]
Con_hub_RBPJmns=Con_hub[2:602,14:26]
Con_hub_RBPJpls=Con_hub[2:602,28:40]

del Con_hub

nuc_DAPT=nuc[2:602,0:20]
nuc_RBPJmns=nuc[2:602,21:41]
nuc_RBPJpls=nuc[2:602,42:62]
nuc_NICD=nuc[2:602,63:83]
del nuc

hub_DAPT=hub[2:602,0:20]
hub_RBPJmns=hub[2:602,21:41]
hub_RBPJpls=hub[2:602,42:62]
del hub

#####

Con_nuc_DAPT=Con_nuc_DAPT.astype(float)
Con_nuc_RBPJmns=Con_nuc_RBPJmns.astype(float)
Con_nuc_RBPJpls=Con_nuc_RBPJpls.astype(float)
Con_nuc_NICD=Con_nuc_NICD.astype(float)

Con_hub_DAPT=Con_hub_DAPT.astype(float)
Con_hub_RBPJmns=Con_hub_RBPJmns.astype(float)
Con_hub_RBPJpls=Con_hub_RBPJpls.astype(float)


nuc_DAPT=nuc_DAPT.astype(float)
nuc_RBPJmns=nuc_RBPJmns.astype(float)
nuc_RBPJpls=nuc_RBPJpls.astype(float)
nuc_NICD=nuc_NICD.astype(float)

hub_DAPT=hub_DAPT.astype(float)
hub_RBPJmns=hub_RBPJmns.astype(float)
hub_RBPJpls=hub_RBPJpls.astype(float)

####

frst_Con_nuc_DAPT=Con_nuc_DAPT[0,:]
frst_Con_nuc_RBPJmns=Con_nuc_RBPJmns[0,:]
frst_Con_nuc_RBPJpls=Con_nuc_RBPJpls[0,:]
frst_Con_nuc_NICD=Con_nuc_NICD[0,:]

frst_Con_hub_DAPT=Con_hub_DAPT[0,:]
frst_Con_hub_RBPJmns=Con_hub_RBPJmns[0,:]
frst_Con_hub_RBPJpls=Con_hub_RBPJpls[0,:]


frst_nuc_DAPT=nuc_DAPT[0,:]
frst_nuc_RBPJmns=nuc_RBPJmns[0,:]
frst_nuc_RBPJpls=nuc_RBPJpls[0,:]
frst_nuc_NICD=nuc_NICD[0,:]

frst_hub_DAPT=hub_DAPT[0,:]
frst_hub_RBPJmns=hub_RBPJmns[0,:]
frst_hub_RBPJpls=hub_RBPJpls[0,:]

###### normalization

Norm_Con_nuc_DAPT=Con_nuc_DAPT[1:600,:]/Con_nuc_DAPT[0,:]
Norm_Con_nuc_RBPJmns=Con_nuc_RBPJmns[1:600,:]/Con_nuc_RBPJmns[0,:]
Norm_Con_nuc_RBPJpls=Con_nuc_RBPJpls[1:600,:]/Con_nuc_RBPJpls[0,:]
Norm_Con_nuc_NICD=Con_nuc_NICD[1:600,:]/Con_nuc_NICD[0,:]

Norm_Con_hub_DAPT=Con_hub_DAPT[1:600,:]/Con_hub_DAPT[0,:]
Norm_Con_hub_RBPJmns=Con_hub_RBPJmns[1:600,:]/Con_hub_RBPJmns[0,:]
Norm_Con_hub_RBPJpls=Con_hub_RBPJpls[1:600,:]/Con_hub_RBPJpls[0,:]

Norm_nuc_DAPT=(nuc_DAPT[1:600,:]-nuc_DAPT[1,:])/nuc_DAPT[0,:]
Norm_nuc_RBPJmns=(nuc_RBPJmns[1:600,:]-nuc_RBPJmns[1,:])/nuc_RBPJmns[0,:]
Norm_nuc_RBPJpls=(nuc_RBPJpls[1:600,:]-nuc_RBPJpls[1,:])/nuc_RBPJpls[0,:]
Norm_nuc_NICD=(nuc_NICD[1:600,:]-nuc_NICD[1,:])/nuc_NICD[0,:]

Norm_hub_DAPT=(hub_DAPT[1:600,:]-hub_DAPT[1,:])/hub_DAPT[0,:]
Norm_hub_RBPJmns=(hub_RBPJmns[1:600,:]-hub_RBPJmns[1,:])/hub_RBPJmns[0,:]
Norm_hub_RBPJpls=(hub_RBPJpls[1:600,:]-hub_RBPJpls[1,:])/hub_RBPJpls[0,:]


#### 





vec = np.linspace(0, 30-0.1, num=599)
vec2 = np.linspace(0, 39.6-0.132, num=599)

ResultMat=np.zeros((4,2,20,3))
HalfLifeMat=np.zeros((4,2,20))

bound1=[0,0,0],[30,1,30]

weights= np.ones_like(vec)
weights[300:599]=100

target_value=0.5

####
for i in range(20):
    print(i)
    print(9<i<15)
    
    if 9<i<15:
        params, covariance = curve_fit(exp_recovary, vec2, Norm_nuc_DAPT[:,i],bounds=bound1, method='trf',sigma=weights, ftol=1e-10, xtol=1e-10, maxfev=10000)
        ResultMat[0,0,i,:]= params
        resultTemp=minimize(lambda x: T_half(x,*params, target_value), 2)
        HalfLifeMat[0,0,i] = resultTemp.x[0]

        params, covariance = curve_fit(exp_recovary, vec2, Norm_nuc_RBPJmns[:,i],bounds=bound1,  method='trf',sigma=weights, ftol=1e-10, xtol=1e-10, maxfev=10000)
        ResultMat[1,0,i,:]= params
        resultTemp=minimize(lambda x: T_half(x,*params, target_value), 2)
        HalfLifeMat[1,0,i] = resultTemp.x[0]
        
        params, covariance = curve_fit(exp_recovary, vec2, Norm_nuc_RBPJpls[:,i],bounds=bound1,  method='trf',sigma=weights, ftol=1e-10, xtol=1e-10, maxfev=10000)
        ResultMat[2,0,i,:]= params
        resultTemp=minimize(lambda x: T_half(x,*params, target_value), 2)
        HalfLifeMat[2,0,i] = resultTemp.x[0]
        
        params, covariance = curve_fit(exp_recovary, vec2, Norm_nuc_NICD[:,i],bounds=bound1,  method='trf',sigma=weights, ftol=1e-10, xtol=1e-10, maxfev=10000)
        ResultMat[3,0,i,:]= params
        resultTemp=minimize(lambda x: T_half(x,*params, target_value), 2)
        HalfLifeMat[3,0,i] = resultTemp.x[0]
        
        params, covariance = curve_fit(exp_recovary, vec2, Norm_hub_DAPT[:,i],bounds=bound1,  method='trf',sigma=weights, ftol=1e-10, xtol=1e-10, maxfev=10000)
        ResultMat[0,1,i,:]= params
        resultTemp=minimize(lambda x: T_half(x,*params, target_value), 2)
        HalfLifeMat[0,1,i] = resultTemp.x[0]
        
        params, covariance = curve_fit(exp_recovary, vec2, Norm_hub_RBPJmns[:,i],bounds=bound1,  method='trf',sigma=weights, ftol=1e-10, xtol=1e-10, maxfev=10000)
        ResultMat[1,1,i,:]= params
        resultTemp=minimize(lambda x: T_half(x,*params, target_value), 2)
        HalfLifeMat[1,1,i] = resultTemp.x[0]
        
        params, covariance = curve_fit(exp_recovary, vec2, Norm_hub_RBPJpls[:,i],bounds=bound1,  method='trf',sigma=weights, ftol=1e-10, xtol=1e-10, maxfev=10000)
        ResultMat[2,1,i,:]= params
        resultTemp=minimize(lambda x: T_half(x,*params, target_value), 2)
        HalfLifeMat[2,1,i] = resultTemp.x[0]
        
    else:
         params, covariance = curve_fit(exp_recovary, vec, Norm_nuc_DAPT[:,i],bounds=bound1, method='trf',sigma=weights, ftol=1e-10, xtol=1e-10, maxfev=10000)
         ResultMat[0,0,i,:]= params
         resultTemp=minimize(lambda x: T_half(x,*params, target_value), 2)
         HalfLifeMat[0,0,i] = resultTemp.x[0]
    
         params, covariance = curve_fit(exp_recovary, vec, Norm_nuc_RBPJmns[:,i],bounds=bound1,  method='trf',sigma=weights, ftol=1e-10, xtol=1e-10, maxfev=10000)
         ResultMat[1,0,i,:]= params
         resultTemp=minimize(lambda x: T_half(x,*params, target_value), 2)
         HalfLifeMat[1,0,i] = resultTemp.x[0]
         
         params, covariance = curve_fit(exp_recovary, vec, Norm_nuc_RBPJpls[:,i],bounds=bound1,  method='trf',sigma=weights, ftol=1e-10, xtol=1e-10, maxfev=10000)
         ResultMat[2,0,i,:]= params
         resultTemp=minimize(lambda x: T_half(x,*params, target_value), 2)
         HalfLifeMat[2,0,i] = resultTemp.x[0]
         
         params, covariance = curve_fit(exp_recovary, vec, Norm_nuc_NICD[:,i],bounds=bound1,  method='trf',sigma=weights, ftol=1e-10, xtol=1e-10, maxfev=10000)
         ResultMat[3,0,i,:]= params
         resultTemp=minimize(lambda x: T_half(x,*params, target_value), 2)
         HalfLifeMat[3,0,i] = resultTemp.x[0]
         
         
         params, covariance = curve_fit(exp_recovary, vec, Norm_hub_DAPT[:,i],bounds=bound1,  method='trf',sigma=weights, ftol=1e-10, xtol=1e-10, maxfev=10000)
         ResultMat[0,1,i,:]= params
         resultTemp=minimize(lambda x: T_half(x,*params, target_value), 2)
         HalfLifeMat[0,1,i] = resultTemp.x[0]
         
         params, covariance = curve_fit(exp_recovary, vec, Norm_hub_RBPJmns[:,i],bounds=bound1,  method='trf',sigma=weights, ftol=1e-10, xtol=1e-10, maxfev=10000)
         ResultMat[1,1,i,:]= params
         resultTemp=minimize(lambda x: T_half(x,*params, target_value), 2)
         HalfLifeMat[1,1,i] = resultTemp.x[0]
         
         params, covariance = curve_fit(exp_recovary, vec, Norm_hub_RBPJpls[:,i],bounds=bound1,  method='trf',sigma=weights, ftol=1e-10, xtol=1e-10, maxfev=10000)
         ResultMat[2,1,i,:]= params
         resultTemp=minimize(lambda x: T_half(x,*params, target_value), 2)
         HalfLifeMat[2,1,i] = resultTemp.x[0]
    
 
mean_NucRBPJmns=np.mean(Norm_nuc_RBPJmns,1)
mean_hubRBPJmns=np.mean(Norm_hub_RBPJmns,1)

mean_NucRBPJpls=np.mean(Norm_nuc_RBPJpls,1)
mean_hubRBPJpls=np.mean(Norm_hub_RBPJpls,1)

mean_NucDAPT=np.mean(Norm_nuc_DAPT,1)
mean_hubDAPT=np.mean(Norm_hub_DAPT,1)

dev_NucRBPJmns=np.std(Norm_nuc_RBPJmns,1)
dev_NucRBPJmns=dev_NucRBPJmns/(20**0.5)

dev_HubRBPJmns=np.std(Norm_hub_RBPJmns,1)
dev_HubsRBPJmns=dev_HubRBPJmns/(20**0.5)

# mobileNucRBPJmns=frst_nuc_RBPJmns-np.max(nuc_RBPJmns[1:-1,:])

sample=19
colorvec=['red','blue','magenta','green']

plt.rcParams.update({'font.size': 14})
plt.plot(vec,Norm_hub_DAPT[:,sample],linestyle='dashed')
y_pred = exp_recovary( vec, *ResultMat[0,1,sample,:])
plt.plot(vec,y_pred,linestyle='dashed', label=HalfLifeMat[0,1,sample],color=colorvec[0])
plt.xlim((-1,30))
plt.show()

plt.plot(vec,mean_hubRBPJpls,linestyle='dashed',color=colorvec[0])
y_pred = exp_recovary( vec,  *ResultMat[3,0,sample,:])
plt.plot(vec,(vec*0 +ResultMat[3,0,sample,1]/2),linestyle='dashed', color=colorvec[2],label=HalfLifeMat[3,0,sample])
plt.plot(vec,y_pred,linestyle='dashed', color=colorvec[2],label=ResultMat[3,0,sample,1]/2)
plt.legend(loc='lower left')
plt.xlim((0,30))

plt.show()


plt.rcParams.update({'font.size': 14})
plt.plot(vec,mean_NucRBPJpls,linestyle='dashed')
params, covariance = curve_fit(exp_recovary, vec, mean_NucRBPJmns,bounds=bound1,  method='trf',sigma=weights, ftol=1e-10, xtol=1e-10, maxfev=10000)
y_pred = exp_recovary( vec, *params)
plt.plot(vec,y_pred,linestyle='dashed', label=HalfLifeMat[0,1,sample],color=colorvec[0])
plt.fill_between(vec, mean_NucRBPJmns-dev_NucRBPJmns, mean_NucRBPJmns+dev_NucRBPJmns,alpha=0.3, color=colorvec[1])
plt.xlim((-1,30))
plt.show()


#plot two averages RBPJ mns Hubs Vs nuc
plt.fill_between(vec, mean_hubRBPJmns-dev_HubsRBPJmns, mean_hubRBPJmns+dev_HubsRBPJmns, color='#ED7D31',label='Hubs')
plt.fill_between(vec, mean_NucRBPJmns-dev_NucRBPJmns, mean_NucRBPJmns+dev_NucRBPJmns, color='#2E75B6', label='Nucleoplasm')
plt.legend(loc='lower right')
plt.xlim((-1,30))
plt.show()

#  for i in range(10):
#      plt.plot(vec,Norm_Con_nuc_DAPT[:,i],linestyle='dashed',color=colorvec[0])
# plt.ylim((0,1.1))
# plt.show()



# for i in [1,3]:
#     y_pred = exp_decay( vec,  C1mean[i],C2mean[i])
#     plt.plot(vec,y_pred,linestyle='dashed',color=colorvec[i])
#     plt.plot(vec,normalized_matrix[ :, i],linestyle='dashed',color=colorvec[i],label=samplevec[i])
#     plt.fill_between(vec, normalized_matrix[ :, i]-ErrorMat_30min[ :, i], normalized_matrix[ :, i]+ErrorMat_30min[ :, i],alpha=0.3,label='_nolegend_',color=colorvec[i])
#     plt.legend(loc='lower left')
# plt.ylim((0,1))
# plt.show()