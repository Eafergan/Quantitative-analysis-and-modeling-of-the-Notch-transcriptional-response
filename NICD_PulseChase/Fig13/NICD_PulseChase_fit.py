from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import multiprocessing
import time


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
    vecC2=np.linspace(0, 0.75, num=51) #background
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
    global normalized_matrix, ErrorMat_30min
    ReturnResult=np.zeros((2,4))
    newMat=newMatGen(normalized_matrix,ErrorMat_30min)
    ReturnResult[:,0]= [0,0]
    ReturnResult[:,1]= C1C2Finder(newMat[:,1])
    ReturnResult[:,2]= [0,0]
    ReturnResult[:,3]= C1C2Finder(newMat[:,3])
    return ReturnResult

df2 = pd.read_excel(r'E:\AOS\1st_daptWashout\OverAll1stNotchDAPTwo.xlsx',sheet_name='Sheet1')
df2AsArray=df2.to_numpy()
flippedArray=np.swapaxes(df2AsArray,0,1)
NumOfTReal=193

samples=['Dll1-FC+, DAPT+','Dll1-FC+','Dll1-FC-','Dll1-FC+, Senexin+']
notch_matrix=np.zeros((NumOfTReal,4 )) 
ErrorMat=np.zeros((NumOfTReal,4 )) 

notch_matrix[:,:]=flippedArray[1:194,24:28]
ErrorMat[:,:]=flippedArray[1:194,30:34]


notch_matrix_30min=np.zeros((186,4))
notch_matrix_30min=notch_matrix[6:194,:]

first_time_points = notch_matrix_30min[0,:]
normalized_matrix = notch_matrix_30min / first_time_points[np.newaxis, :]
ErrorMat_30min=np.zeros((186,4))
ErrorMat_30min=ErrorMat[6:194,:]/ first_time_points[np.newaxis, :]


vec = np.linspace(0, (NumOfTReal-1)*5-30, num=NumOfTReal-6)
fullVec=np.linspace(0, (NumOfTReal-1)*5, num=NumOfTReal)


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
    results[i,:]=IterationFunction(i)
    if i%10==0:
        print(i)
    
C1mean=np.zeros((4)) 
C1std=np.zeros((4)) #0 is mean 1 is std
C2mean=np.zeros((4)) #0 is mean 1 is std
C2std=np.zeros((4))
Halflives_lowest95=np.zeros((4)) 
Halflives_Highest95=np.zeros((4)) 

HalfLives=(np.log(2)/(results[:,0,:]))


for i in range(4):  
    C1mean[i]=np.mean(results[:,0,i])
    C1std[i]=np.std(results[:,0,i])
    C2mean[i]=np.mean(results[:,1,i])
    C2std[i]=np.std(results[:,1,i])
    Halflives_lowest95[i]=(np.percentile(HalfLives[:,i],0.5))
    Halflives_Highest95[i]=(np.percentile(HalfLives[:,i],99.5))

HalfLifeAVG95=(Halflives_lowest95+Halflives_Highest95)*0.5
HalfLifeDEV95=(Halflives_lowest95-Halflives_Highest95)*0.5





HalfLives=(np.log(2)/(results[:,0,:]))

colorvec=['red','blue','magenta','green']
samplevec=['Dll1-Fc+, DAPT+','Dll1-Fc+','Dll1-Fc-','Dll1-Fc+, SenexinA+']
plt.rcParams.update({'font.size': 14})
for i in [2,1,0,3]:
    if i==1 or i==3:
        y_pred = exp_decay( vec,  C1mean[i],C2mean[i])
        plt.plot(vec+30,y_pred*first_time_points[i],linestyle='dashed',color=colorvec[i])
    #plt.plot(fullVec,notch_matrix[ :, i],linestyle='dashed',color=colorvec[i],label=samplevec[i])
    plt.fill_between(fullVec, notch_matrix[ :, i]-ErrorMat[ :, i], notch_matrix[ :, i]+ErrorMat[ :, i],alpha=0.3,label=samplevec[i],color=colorvec[i])
    plt.legend(loc='upper right')
plt.xlim((0,930))
plt.ylim((0,5))
plt.show()


plt.hist(HalfLives[:,1], bins=30, edgecolor='black')  # Adjust the number of bins as needed
plt.grid(True)
plt.show()


# for i in [1,3]:
#     y_pred = exp_decay( vec,  C1mean[i],C2mean[i])
#     plt.plot(vec,y_pred,linestyle='dashed',color=colorvec[i])
#     plt.plot(vec,normalized_matrix[ :, i],linestyle='dashed',color=colorvec[i],label=samplevec[i])
#     plt.fill_between(vec, normalized_matrix[ :, i]-ErrorMat_30min[ :, i], normalized_matrix[ :, i]+ErrorMat_30min[ :, i],alpha=0.3,label='_nolegend_',color=colorvec[i])
#     plt.legend(loc='lower left')
# plt.ylim((0,1))
# plt.show()



# y_pred = exp_decay( vec,  C1mean,C2DAPTmean, C3DAPTmean)
# plt.plot(vec,y_pred,linestyle='dashed',label='+Dll1-Fc, +DAPT')
# plt.fill_between(vec, normalized_matrix[ :, 0]-ErrorMat[ :, 0], normalized_matrix[ :, 0]+ErrorMat[ :, 0],alpha=0.3)
# y_pred = exp_decay( vec,  C1mean,C2Dllmean, C3Dllmean)
# plt.plot(vec,y_pred,linestyle='dashed',label='+Dll1-Fc')
# plt.fill_between(vec, normalized_matrix[ :, 1]-ErrorMat[ :, 1], normalized_matrix[ :, 1]+ErrorMat[ :, 1],alpha=0.3)
# y_pred = exp_decay( vec,  C1mean,C2PBSstd, C3PBSmean)
# plt.plot(vec,y_pred,linestyle='dashed',label='- Dll1-Fc')
# plt.fill_between(vec, normalized_matrix[ :, 2]-ErrorMat[ :, 2], normalized_matrix[ :, 2]+ErrorMat[ :, 2],alpha=0.3)
# plt.ylim((0,1))
# plt.legend()
# plt.show()


plt.plot(vec,exp_decay(vec,0.00545197,0.3),linestyle='dashed',label='fit DAPT')
plt.plot(vec,normalized_matrix[:,1],linestyle='-', color='red', label='DAPT')
plt.legend()
plt.ylim((0,1))
plt.show()

# plt.plot(vec,exp_decay(vec,C1mean,C2Dllmean,C3Dllmean),linestyle='dashed',label='fit Dll')
# plt.plot(vec,normalized_matrix[:,1],linestyle='-', color='red', label='Dll')
# plt.legend()
# plt.ylim((0,1))
# plt.show()

# plt.plot(vec,exp_decay(vec,C1mean,C2PBSmean,C3PBSmean),linestyle='dashed',label='fit PBS')
# plt.plot(vec,normalized_matrix[:,2],linestyle='-', color='red', label='PBS')
# plt.legend()
# plt.ylim((0,1))
# plt.show()'

# Example data (replace this with your actual data)
