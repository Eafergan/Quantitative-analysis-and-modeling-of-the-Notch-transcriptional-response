# Fit fluorescence against known JFX554 concentrations.
%reset -f


from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def Linear_Regression(x, C1,C2): # C1 is the slope; C2 is the intercept.
    return  C1*x + C2


def newVecGen(means,StErr):
    newVec=np.zeros_like(means)
    # Draw one fluorescence value for each of the four concentrations.
    for i in range(4):
        newVec[i]=np.random.normal(means[i], StErr[i])
    return newVec



# Dye concentrations are in nM; means and StErr are fluorescence measurements.
vec=[600,300,150,75]
# Store slope and intercept in rows, with one column per fit.
ReturnResult=np.zeros((2,2000))
means=[4.0115,2.2875,1.57625,1.26525]
StErr=[0.082357832,	0.035306692,	0.025967227,	0.012315514]


# Fit 2,000 random datasets drawn from the measured means and errors.
for i in range(2000):
    newVec=newVecGen(means,StErr)
    params, covariance = curve_fit(Linear_Regression, vec, newVec,bounds=([0,0],[0.1,2]))
    ReturnResult[:,i]= params
    del newVec,params



# Keep the 2.5th and 97.5th percentiles for the slope and intercept.
Mean_lowestC1=(np.percentile(ReturnResult[0,:],2.5))
Mean_HighestC1=(np.percentile(ReturnResult[0,:],97.5))

Mean_lowestC2=(np.percentile(ReturnResult[1,:],2.5))
Mean_HighestC2=(np.percentile(ReturnResult[1,:],97.5))
