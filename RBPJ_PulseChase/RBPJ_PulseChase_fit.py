from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import multiprocessing
import time

#values from previus fit are DAPT[0.00033995,0.0242347,0.885518], DLL:[0.0004781,0.02715,0.06908] PBS:[0.00030004,0.06908742,0.88663632]

def exp_decay(x, C1,C2): #C1 is ratio, C2 is Lambda for 1 and C3 is for 2, different than in the PPT
    return  (1-C2)*np.exp(-x*C1) + C2



def scorer(C1,C2,y):
    global vec
    distance=0
    y_prediction=exp_decay(vec,C1,C2)
    distance_vec=y_prediction-y
    distance=sum(x ** 2 for x in (distance_vec))
    return distance

def C1C2Finder(yvec):
    global vec
    vecC1=np.linspace(0, 1, num=11)
    vecC2=np.linspace(0, 0.75, num=51) #backgroundw
    resultsC2=np.zeros((51))
    resultsC1score=np.zeros((11))
    resultsC1index=np.zeros((16))
    for i in range(len(resultsC1index)):
        for c1 in range(11):
            for c2 in range(51):
                resultsC2[c2]=scorer(vecC1[c1],vecC2[c2], yvec)
            resultsC1score[c1]=scorer(vecC1[c1],vecC2[np.argmin(resultsC2)] , yvec)
        resultsC1index[i]=np.argmin(resultsC1score)
        if resultsC1index[i]>8:
            upperC1=vecC1[-1]
            lowerC1=vecC1[5]
            vecC1=np.linspace(lowerC1, upperC1, num=11)
        else:
            if resultsC1index[i]<3:
                upperC1=vecC1[5]
                lowerC1=vecC1[0]
                vecC1=np.linspace(lowerC1, upperC1, num=11)
            else:
                upperC1=vecC1[int(resultsC1index[i])+2]
                lowerC1=vecC1[int(resultsC1index[i])-3]
                vecC1=np.linspace(lowerC1, upperC1, num=11)
    return vecC1[np.argmin(resultsC1score)],vecC2[np.argmin(resultsC2)]

 
def newMatGen(normalized_matrix,ErrorMat):
    newMat=np.zeros_like(normalized_matrix)
    for sample in range(len(normalized_matrix[0,:])):
        for i in range(len(normalized_matrix[:,0])):
            newMat[i,sample]=np.random.normal(normalized_matrix[i,sample], ErrorMat[i,sample])
    return newMat

def IterationFunction(i):
    global normalized_matrix, ErrorMat
    ReturnResult=np.zeros((2,4))
    newMat=newMatGen(normalized_matrix,ErrorMat)
    ReturnResult[:,0]= C1C2Finder(newMat[:,0])
    ReturnResult[:,1]= C1C2Finder(newMat[:,1])
    ReturnResult[:,2]= C1C2Finder(newMat[:,2])
    ReturnResult[:,3]= C1C2Finder(newMat[:,3])
    return ReturnResult

df2 = pd.read_excel(r'D:\RBPJ\2nd_RBPJ_DAPTwo_48hrs\2nd_RBPJ_DAPTwo_48hrs_OverAll.xlsx',sheet_name='Sheet3_Ch2Norm')
df2AsArray=df2.to_numpy()
flippedArray=np.swapaxes(df2AsArray,0,1)
NumOfTReal=145

normalized_matrix=np.zeros((NumOfTReal,4 )) 
ErrorMat=np.zeros((NumOfTReal,4 )) 

normalized_matrix[:,:]=flippedArray[1:146,24:28]
ErrorMat[:,:]=flippedArray[1:146,30:34]

vec = np.linspace(0, (NumOfTReal-1)*20, num=NumOfTReal)



# start_time_par= time.time() #took 2 min to run 5 iterations 5.5 min for 20 iterations, 110 min to run 500
# if __name__ == "__main__":
#     num_cores = multiprocessing.cpu_count()
#     pool = multiprocessing.Pool(num_cores)

#     num_iterations = 500  
#     results_par = pool.map(IterationFunction, range(num_iterations))
    
#     pool.close()
#     pool.join()

# end_time_par= time.time()
# runtime=end_time_par-start_time_par


num_iterations=500
results=np.zeros((num_iterations,2,4))



for i in range(num_iterations):
    results[i,:,:]=IterationFunction(i)
    if i%10==0:
        print(i)
    
C1mean=np.zeros((4)) 
C1std=np.zeros((4)) #0 is mean 1 is std
C2mean=np.zeros((4)) #0 is mean 1 is std
C2std=np.zeros((4))

for i in range(4):  
    C1mean[i]=np.mean(results[:,0,i])
    C1std[i]=np.std(results[:,0,i])
    C2mean[i]=np.mean(results[:,1,i])
    C2std[i]=np.std(results[:,1,i])
 


HalfLifeplus=(np.log(2)/(C1mean-C1std))
HalfLifeminus=(np.log(2)/(C1mean+C1std))
HalfLifeAVG=(HalfLifeplus+HalfLifeminus)*0.5
HalfLifeDEV=(HalfLifeplus-HalfLifeminus)*0.5


colorvec=['red','blue','magenta','green']
samplevec=['Dll1-FC+, DAPT+','Dll1-FC+','Dll1-FC-','Dll1-FC+, Senexin+']

plt.rcParams.update({'font.size': 14})
y_pred = exp_decay( vec,  C1mean[2],C2mean[2])
plt.plot(vec,y_pred,linestyle='dashed',label=samplevec[2], color=colorvec[2])
plt.fill_between(vec, normalized_matrix[ :, 2]-ErrorMat[ :, 2], normalized_matrix[ :, 2]+ErrorMat[ :, 2],alpha=0.3, color=colorvec[2])
y_pred = exp_decay( vec,  C1mean[1],C2mean[1])
plt.plot(vec,y_pred,linestyle='dashed',label=samplevec[1], color=colorvec[1])
plt.fill_between(vec, normalized_matrix[ :, 1]-ErrorMat[ :, 1], normalized_matrix[ :, 1]+ErrorMat[ :, 1],alpha=0.3, color=colorvec[1])
y_pred = exp_decay( vec,  C1mean[0],C2mean[0])
plt.plot(vec,y_pred,linestyle='dashed',label=samplevec[0],color=colorvec[0])
plt.fill_between(vec, normalized_matrix[ :, 0]-ErrorMat[ :, 0], normalized_matrix[ :, 0]+ErrorMat[ :, 0],alpha=0.3, color=colorvec[0])
y_pred = exp_decay( vec,  C1mean[3],C2mean[3])
plt.plot(vec,y_pred,linestyle='dashed',label=samplevec[3], color=colorvec[3])
plt.fill_between(vec, normalized_matrix[ :, 3]-ErrorMat[ :, 3], normalized_matrix[ :, 3]+ErrorMat[ :, 3],alpha=0.3, color=colorvec[3])
plt.ylim((0,1))
plt.xlim((0,2880))
plt.legend(loc=1,prop={'size': 12})
plt.show()

plt.rcParams.update({'font.size': 14})
y_pred = exp_decay( vec,  C1mean[2],C2mean[2])
plt.plot(vec,y_pred,linestyle='dashed',label=samplevec[2], color=colorvec[2])
plt.fill_between(vec, normalized_matrix[ :, 2]-ErrorMat[ :, 2], normalized_matrix[ :, 2]+ErrorMat[ :, 2],alpha=0.3, color=colorvec[2])
y_pred = exp_decay( vec,  C1mean[1],C2mean[1])
plt.plot(vec,y_pred,linestyle='dashed',label=samplevec[1], color=colorvec[1])
plt.fill_between(vec, normalized_matrix[ :, 1]-ErrorMat[ :, 1], normalized_matrix[ :, 1]+ErrorMat[ :, 1],alpha=0.3, color=colorvec[1])
plt.ylim((0.2,1))
plt.xlim((0,2000))
plt.legend(loc=1,prop={'size': 12})
plt.show()

# plt.plot(vec,exp_decay(vec,C1mean,C2DAPTmean,C3DAPTmean),linestyle='dashed',label='fit DAPT')
# plt.plot(vec,normalized_matrix[:,0],linestyle='-', color='red', label='DAPT')
# plt.legend()
# plt.ylim((0,1))
# plt.show()

# plt.plot(vec,exp_decay(vec,C1mean,C2Dllmean,C3Dllmean),linestyle='dashed',label='fit Dll')
# plt.plot(vec,normalized_matrix[:,1],linestyle='-', color='red', label='Dll')
# plt.legend()
# plt.ylim((0,1))
# plt.show()

# plt.plot(vec,exp_decay(vec,C1mean,C2PBSmean,C3PBSmean),linestyle='dashed',label='fit PBS')
# plt.plot(vec,normalized_matrix[:,2],linestyle='-', color='red', label='PBS')
# plt.legend()
# plt.ylim((0,1))
# plt.show()