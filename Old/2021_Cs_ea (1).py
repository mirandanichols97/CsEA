# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 00:51:30 2021

@author: J Edu
"""

import matplotlib.pyplot as plt
# import csv
import numpy as np
# pylab.rcParams['figure.figsize'] = (13, 8)

import pandas as pd

from scipy.optimize import curve_fit

def Wigner(x,a,b,ea):
    return a + b*np.sqrt(x-ea+np.abs(x-ea))
    
exp_path = 'Users/xnimir/Desktop/'
file_name = 'counter_prop_WignerS'
save_img = False

# exp_sheets = ['Counter Feb 11_new'] # put here all co and counter prop measurements, (counters are even, starting from 0)
exp_sheets = ['Counter Feb 13'] # put here all co and counter prop measurements, (counters are even, starting from 0)
columns = ['Signal',' WL red (nm)','Ion current (pA)','Laser power red (mW)','std dev','pow dev']
n = 3
c= 299792458
hc = 1239.84193 #eV/nm, this is the factor to convert nm -> eV
ea_s = []
err_ea_s = []

def np_array(lista):
    li0 = lista.to_numpy()
    li = li0.flatten()
    return np.array(li[:])

for k,sheet in enumerate(exp_sheets):
    d_params = pd.read_excel(exp_path + 'Cs_EA_new.xlsx'
                             ,sheet_name = sheet) # need to fill in the excel sheet
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
    
    if k == 0:
        x1 = ens[5:-2]
        y1 = s_n[5:-2]
        err_y = var_s[5:-2]
    else:
        x1 = ens[:]
        y1 = s_n[:]
        err_y = var_s[:]        
    
    if divmod(k,2)[1] == 1: 
        pars = [1.7,300,1.9268]
        xarray=np.arange(pars[2]-400e-6,pars[2]+200e-6,1e-7)
        doppler_est = -600*1e-6 
    elif divmod(k,2)[1] == 0:
        pars = [1,100,1.92563]
        xarray=np.arange(pars[2]-400e-6,pars[2]+200e-6,1e-7)
        # print (xarray)
        doppler_est = 600*1e-6

    setColor1='black'
    setColor2=[1,0,0]
    # xarray=np.arange(x1[0],x1[len(x1)-1]+0.0005,0.0000001)
    y_func = Wigner(xarray,*pars)
    # bound_par = ([0,0,ea-0.01],[1,200,ea+0.01])
    #pl2, covl2 = curve_fit(Wigner, x1, y1,[a,b,ea],bounds=bound_par) # for bound fit
    try:
        pl0, covl0 = curve_fit(Wigner, x1, y1,[*pars],sigma = err_y, absolute_sigma=True,maxfev=1000)
    # print (pl0)
        plt.figure(k)
        plt.errorbar(x1, y1, yerr=err_y, fmt='.k',color='blue', capthick=0.5,capsize=5,elinewidth=0.5)
        plt.plot(xarray,Wigner(xarray, *pl0),label='fit',color=setColor1)
        plt.plot(xarray,y_func, linestyle = '--',alpha = 0.6)
        plt.scatter(x1,y1,s=7)
        plt.xlabel("Photon energy (eV)")
        plt.ylabel("counts(normalized)")
        st_error = np.sqrt(np.diag(covl0))[2]
        err_ea_s.append(st_error)
        plt.title(sheet)
        middle_trans = 1.454620680
        Doppler = 600*1e-6
        electron_affinity = abs(pl0[2]-middle_trans)
        wl_ea = hc/pl0[2]
        wls = hc/x1
        ea_s.append(electron_affinity)
        print (round(electron_affinity + doppler_est,7), 'eV', 'statistichal error:', str(round(st_error,7)))
        if save_img:
            plt.savefig(exp_path + file_name + '_{}.png'.format(sheet), dpi = 300, bbox_inches='tight')
        plt.show()
    except Exception as e:
        print (e)
    plt.figure(k+20)
    plt.scatter(wls,y1, s=5)
    plt.errorbar(wls, y1, yerr=err_y, fmt='.k',color='blue', capthick=0.5,capsize=5,elinewidth=0.5)
    plt.axvline(x = wl_ea, color = 'pink', linestyle = '--')
    plt.xlabel("wavelength (nm)")
    plt.ylabel("counts(normalized)")
    plt.xlim([wl_ea-0.03,wl_ea+0.03])
    plt.text(wls[2], y1[3], exp_sheets[k] + ','+ str(round(wl_ea,4)) + ' nm', bbox=dict(facecolor='red', alpha=0.3))

    if save_img:
        plt.savefig(exp_path + file_name + 'wl_{}.png'.format(sheet), dpi = 300, bbox_inches='tight')
    plt.show()


#%%


#%%
################################### fitting part

file_name = 'Cs-_ea_s_new'
save_img = True

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
    
x2.append(2)
x2.append(3)
y2.append(471611.5)
y2.append(471626)
err2.append(1)
err2.append(25)
plt.figure(10)
plt.scatter(x2,y2, s=7)
plt.errorbar(x2,y2,yerr = err2,fmt='.k',color='blue', capthick=0.5,capsize=5,elinewidth=0.5)
plt.grid(lw = 0.3)
plt.ylim([min(y2)-15e0,max(y2)+26e0])
plt.text(x2[0], y2[0]+3, exp_sheets[0] + ','+ exp_sheets[1], bbox=dict(facecolor='red', alpha=0.3))
plt.text(x2[1]+0.1, y2[1]-5, exp_sheets[1] + ','+ exp_sheets[2], bbox=dict(facecolor='red', alpha=0.3))
plt.text(x2[2]-0.4, y2[2], 'Thesis', bbox=dict(facecolor='red', alpha=0.3))
plt.text(x2[3]-0.9, y2[3], 'Anderson et. al.', bbox=dict(facecolor='red', alpha=0.3))
plt.title('Geo. mean')
plt.ylabel('energy (micro eV)')
plt.xlabel('no. experiment')
if save_img:
    plt.savefig(exp_path + file_name + '.png', dpi = 300, bbox_inches='tight')
plt.show()

x3 = np.arange(0,2,1)
y3 = ea_s[::2]
err3 = err_ea_s[::2]
y3 = np.array(y3)
err3 = np.array(err3)
plt.figure(30)
plt.scatter(x3,y3,s=5)
plt.errorbar(x3, y3, yerr=err3, fmt='.k',color='blue', capthick=0.5,capsize=5,elinewidth=0.5)
plt.ylim([y3[0]-0.00002,y3[0]+0.000005])
plt.text(x3[0]+0.1, y3[0], exp_sheets[0], bbox=dict(facecolor='red', alpha=0.3))
plt.text(x3[1]-0.3, y3[1], exp_sheets[2], bbox=dict(facecolor='red', alpha=0.3))
plt.title('EA for counter')
plt.ylabel('energy (eV)')
plt.xlabel('no. experiment')
if save_img:
    plt.savefig(exp_path + file_name + '_counter.png', dpi = 300, bbox_inches='tight')
plt.show()
