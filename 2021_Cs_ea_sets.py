# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 00:51:30 2021

@author: J Edu Hola! 
"""

import matplotlib.pyplot as plt
import numpy as np
# pylab.rcParams['figure.figsize'] = (13, 8)
import pandas as pd
from scipy.optimize import curve_fit

def Wigner(x,a,b,ea):
    return a + b*np.sqrt(x-ea+np.abs(x-ea))
    
exp_path = '/Users\Edu\Documents\PhD programme\Cs exp GU/'
file_name = 'WignerS' # name for the saved plot, png, dpi = 300
save_img = False
xl_name = 'Cs_EA_set4_Feb19.xlsx'

# put here the xl sheets from the corresponding sets
exp_sheets = ['Co1','Counter1','Co2','Counter2','Co3']
columns = ['Signal',' WL red (nm)','Ion current (pA)','Laser power red (mW)','std dev','pow dev']

c= 299792458
hc = 1239.84193 #eV/nm, this is the factor to convert nm -> eV
ea_s = []
err_ea_s = []

def np_array(lista):
    li0 = lista.to_numpy()
    li = li0.flatten()
    return np.array(li[:])

for k,sheet in enumerate(exp_sheets):
    d_params = pd.read_excel(exp_path + xl_name,sheet_name = sheet)
    data = []
    for col in columns:
        df_s = pd.DataFrame(d_params,columns = [col])     
        data.append(np_array(df_s))
    
    for j,counts in enumerate(data[0]):
        if counts == 0:
            for i in range(len(data)):
                np.delete(data[i],j) # this takes care of any zeros involved
                
    ens = hc/data[1]
    s_n = data[0]/(data[2]*data[3])
    var_s0 = (np.sqrt(data[0])/(data[2]*data[3]))**2 + (data[0]/(data[2]**2*data[3])*data[4])**2 + (data[0]/(data[2]*data[3]**2)*data[5])**2
            
    var_s = np.sqrt(var_s0)
    # d = data[0].sort()
    # print (len(data[1]),d)
    
    if k == 1:
        x1 = ens[:]
        y1 = s_n[:]
        err_y = var_s[:]
    else:
        x1 = ens[:]
        y1 = s_n[:]
        err_y = var_s[:]        

# this defines initial parameters, if Co is the first element, set this value to 0, otherwise to 1
    if divmod(k,2)[1] == 0:   
        pars = [1.7,300,1.9268]
        xarray=np.arange(pars[2]-100e-6,pars[2]+200e-6,1e-7)
        doppler_est = -600*1e-6 
# if Co is the first element, set this value to 1, otherwise to 0
    elif divmod(k,2)[1] == 1:
        pars = [1,100,1.92563]
        xarray=np.arange(pars[2]-100e-6,pars[2]+200e-6,1e-7)
        # print (xarray)
        doppler_est = 600*1e-6

    
    y_func = Wigner(xarray,*pars) # definition of the trial Wigner curve
    bound_par = ([0,0,pars[2]-0.0002],[5,5000,pars[2]+0.0001])
    #pl2, covl2 = curve_fit(Wigner, x1, y1,[a,b,ea],bounds=bound_par) # for bound fit
    try:
        pl0, covl0 = curve_fit(Wigner, x1, y1,[*pars],bounds=bound_par,sigma = err_y, absolute_sigma=True,maxfev=1000)
    # print (pl0)
        plt.figure(k)
        plt.errorbar(x1, y1, yerr=err_y, fmt='.k',color='blue', capthick=0.5,capsize=5,elinewidth=0.5)
        plt.plot(xarray,Wigner(xarray, *pl0),label='fit',color='black')
        plt.plot(xarray,y_func, linestyle = '--',alpha = 0.6) # plotting of the trial Wigner curve
        plt.scatter(x1,y1,s=7)
        plt.xlabel("Photon energy (eV)")
        plt.ylabel("counts(normalized)")
        st_error = np.sqrt(np.diag(covl0))[2] # Takes cov Matrix and gives the Std of the parameters
        err_ea_s.append(st_error)
        plt.title(sheet)
        middle_trans = 1.454620680 # energy diff between Ground state and 6p in eV
        # Doppler = 600*1e-6 # this is just to estimate the Doppler shift if needed
        electron_affinity = abs(pl0[2]-middle_trans)
        ea_s.append(electron_affinity)
        plt.text(x1[0],y1[7],'EA: ' + str(round(electron_affinity,7)) + ' eV', bbox=dict(facecolor='red', alpha=0.3))
        print (round(electron_affinity,7), 'eV', 'statistichal error:', str(round(st_error,7)))
        if save_img:
            plt.savefig(exp_path + file_name + '_{}.png'.format(sheet), dpi = 300, bbox_inches='tight')
        plt.show()
    except Exception as e:
        print (e)
        
# # this part is for plottin the curve with wavelengths
    wls = hc/x1
    wl_ea = hc/pl0[2] # gives the threshold wavelength
    plt.figure(k+20)
    plt.scatter(wls,y1, s=5)
    plt.errorbar(wls, y1, yerr=err_y, fmt='.k',color='blue', capthick=0.5,capsize=5,elinewidth=0.5)
    plt.axvline(x = wl_ea, color = 'pink', linestyle = '--')
    plt.xlabel("wavelength (nm)")
    plt.ylabel("counts(normalized)")
    plt.xlim([wl_ea-0.03,wl_ea+0.07])
    plt.text(wls[2], y1[3], exp_sheets[k] + ','+ str(round(wl_ea,4)) + ' nm', bbox=dict(facecolor='red', alpha=0.3))

    if save_img:
        plt.savefig(exp_path + file_name + 'wl_{}.png'.format(sheet), dpi = 300, bbox_inches='tight')
    plt.show()



#%%
################################### Calc of Geo Mean of EA's

x2 = []
y2 = []
err2 = []
for i in range(len(exp_sheets)-1):
    ele_aff = np.sqrt(ea_s[i]*ea_s[i+1]) 
    exp_error = (1/2)*(ea_s[i]*ea_s[i+1])**(-1/2)*np.sqrt((ea_s[i]*err_ea_s[i+1])**2 + (ea_s[i]*err_ea_s[i+1])**2)
    x2.append(i)
    y2.append(ele_aff*1e6)
    err2.append(exp_error*1e6)
    print ('Cs- electron aff:', str(round(ele_aff,7)), 'error: ', str(round(exp_error,7)))
 

#%%
############### With this part you can see the plot of EA's and also compared them with Literature value

from statsmodels.stats.weightstats import DescrStatsW

save_img = False
file_name = 'Cs-_ea_overDays'

x_01 = np.arange(1,3,1)
EA_1 = np.array([0.4716150,0.4716083])*1e6
EA_err_1 = np.array([0.6,1])
x_02 = np.arange(3,5,1)
EA_2 = np.array([0.4716198,0.4716201])*1e6
EA_err_2 = np.array([0.8,0.7])
x_03 = np.arange(5,7,1)
EA_3 = np.array([0.4716198,0.4716227])*1e6
EA_err_3 = np.array([0.3,1.4])
x_04 = np.arange(7,10,1)
EA_4 = np.array([0.4716217,0.4716233,0.4716227])*1e6
EA_err_4 = np.array([0.36,0.56,0.26])
x_05 = np.arange(10,14,1)
EA_5 = np.array([0.4716183,0.4716189,0.4716212,0.4716210])*1e6
EA_err_5 = np.array([0.9,0.3,0.4,0.8])

EA_T =  EA_2.tolist() + EA_3.tolist()+ EA_4.tolist() + EA_5.tolist()
EA_err_T = EA_err_2.tolist() + EA_err_3.tolist()+ EA_err_4.tolist()+ EA_err_5.tolist()
np.average(EA_T,weights = EA_err_T)
w_stats = DescrStatsW(EA_T, weights=EA_err_T, ddof=0)

EA_6 = np.array([0.4716115,0.471626])*1e6
EA_err_6 = np.array([1,25])
x_06 = np.arange(14,16,1)

plt.figure(40)
# plt.scatter(x_01,EA_1, color = 'blue',label = 'Feb-11, Feb-12')
# plt.errorbar(x_01, EA_1, yerr=EA_err_1, fmt='.k',color='blue', capthick=0.5,capsize=5,elinewidth=0.5)
plt.scatter(x_02,EA_2, color = 'black',label = 'Feb-16')
plt.errorbar(x_02, EA_2, yerr=EA_err_2, fmt='.k',color='black', capthick=0.5,capsize=5,elinewidth=0.5)
plt.scatter(x_03,EA_3, color = 'orange',label = 'Feb-17')
plt.errorbar(x_03, EA_3, yerr=EA_err_3, fmt='.k',color='orange', capthick=0.5,capsize=5,elinewidth=0.5)
plt.scatter(x_04,EA_4, color = 'purple',label = 'Feb-18')
plt.errorbar(x_04, EA_4, yerr=EA_err_4, fmt='.k',color='purple', capthick=0.5,capsize=5,elinewidth=0.5)
plt.scatter(x_05,EA_5, color = 'red',label = 'Feb-19')
plt.errorbar(x_05, EA_5, yerr=EA_err_5, fmt='.k',color='red', capthick=0.5,capsize=5,elinewidth=0.5)

plt.scatter(x_06,EA_6, color = 'green',label = 'Thesis, Andersen et.al.')
plt.errorbar(x_06, EA_6, yerr=EA_err_6, fmt='.k',color='green', capthick=0.5,capsize=5,elinewidth=0.5)
plt.axhline(y = w_stats.mean,color = 'blue', linestyle = '--')
plt.axhline(y =w_stats.mean + w_stats.std,color = 'black', linestyle = '--', alpha = 0.3)
plt.axhline(y =w_stats.mean - w_stats.std,color = 'black', linestyle = '--', alpha = 0.3)
plt.text(x_03[1]+2.5, EA_3[1] + 170*1e-1, 
         'EA: ' + str(round(w_stats.mean*1e-6,7)) + '(' + str(int(w_stats.std*1e1))+')' + ' eV', bbox=dict(facecolor='red', alpha=0.3))
plt.legend(loc = 'best')
plt.grid(lw = 0.3)
plt.title('Geo. mean')
plt.ylabel('energy (micro eV)')
plt.xlabel('no. experiment')
if save_img:
    plt.savefig(exp_path + file_name + '.png', dpi = 300, bbox_inches='tight')
plt.show()