"""
Galah
Ultra 
Enigmatic
Spectroscopic
Script

last updated: 25-10-17

@author: Jane Lin (u5027368@anu.edu.au)

*The RV portion of GUESS, splits CCDs 2 and 3 into 4 parts and calculates RV for each segment, used for testing wavelength calibration issues
*If the cross correlation fails, returns np.nan
*Also returns the sigma of the segments
*Values are now 3 dp 


"""

import pyfits
import scipy.signal as sig
import numpy as np
from scipy.optimize import leastsq
import glob
import pdb
import csv
from itertools import izip
from scipy.interpolate import UnivariateSpline
import heapq
import matplotlib.pyplot as plt
from numpy.linalg import lstsq
import pickle
import os, sys
from scipy.stats import sigmaclip

guess_params=[]
for line in open("guess_param.txt"):
    li=line.strip()
    if not li.startswith("#"):
        guess_params.append(li)
#-------------------------------saving locations-----------------------#

folder=guess_params[0] 
saveloc=guess_params[1] 
folderloc=guess_params[2] 
param_output=guess_params[3] 

wlgrid1=guess_params[4] 
wlgrid2=guess_params[5] 
wlgrid3=guess_params[6] 

rvmodel1=guess_params[7] 
rvmodel2=guess_params[8] 
rvmodel3=guess_params[9] 

plscale1=guess_params[10] 
plscale2=guess_params[11] 
plscale3=guess_params[12] 
plscale4=guess_params[13] 

continuum_reg=guess_params[14] 
param_grid=guess_params[15] 

fit_degree1=int(guess_params[19])
fit_degree2=int(guess_params[20])
fit_degree3=int(guess_params[21])
fit_degree4=int(guess_params[22])

#----------------------------------------------------------------------#

#wavelength grid
wave1=np.loadtxt(wlgrid1)
wave2=np.loadtxt(wlgrid2)
wave3=np.loadtxt(wlgrid3)

if not os.path.exists(folderloc):
    os.makedirs(folderloc)
    print 'created %s' %folderloc

sigma=5 #for the median filter 
telluric_region=int(guess_params[17]) #default=7707A, cutting off anything before this
dir_str=guess_params[18]
#dir_str='/priv/miner3/galah/2dfdr_reductions/aaorun/614'

bad_weights=[]#keeps track of stars with bad cc
sn_low=[]#low sn

def clean(filename):
    
    ''' 
    inputs:
    1. filename of the unnormalised 2dfdr spectrum, be like 
    '/priv/miner3/galah/2dfdr_reductions/aaorun/614/150606/1506060029013971.fits' (ccd1)
    
    outputs:
    1. interpolated, median filtered flux
    2. the common wl grid
    3. s/n
    4. median filtered flux (un-interploated)
    5. original wl grid, without nans                 
    '''
    
    if filename.split('.')[-2].endswith('1'):
        grid=wave1
    if filename.split('.')[-2].endswith('2'):
        grid=wave2
    if filename.split('.')[-2].endswith('3'):
        grid=wave3
    
    #starno=filename.split('_')[-1].split('.')[0]#galahic
    starno=filename.split('.')[:-1]
    hdu=pyfits.open(filename)
    h = pyfits.getheader(filename)
    x = np.arange(h['NAXIS1']) + 1 - h['CRPIX1']
    wave = h['CDELT1']*x + h['CRVAL1']
    d=hdu[0].data  
    #getting rid of the nan values at the beginning of the ccd
    nan_end=np.where(np.isnan(np.array(hdu[0].data)[::-1][:2000]))
    nan_beg=np.where(np.isnan(np.array(hdu[0].data)[:2000]))
    try:
        d=hdu[0].data[:-1*(nan_end[0][-1]+1)]
        wave=wave[:-1*(nan_end[0][-1]+1)]
    except:
        pass
    try:
        d=d[nan_beg[0][-1]+1:]
        wave=wave[nan_beg[0][-1]+1:]
    except:
        pass   
    nans2=np.where(np.isnan(d))[0] #for nans in the middle. but WHY ARE THERE NANS IN THE MIDDLE!?!!? GARH!
    if len(nans2)!=0:
        d[nans2]=np.mean(d[np.where(np.isnan(d)==False)])
    #median filter smoothing and calculating s/n
      
    signal = np.median(d)
    noise = 1.4826/np.sqrt(2)*np.median(np.abs(d[1:] - d[:-1]))
    snr = signal / noise
    med_d = sig.medfilt(d,5) #5 pixels window
    w = np.where(abs(med_d - d) > sigma * noise)[0]
    d[w]=med_d[w] 
    sn=True

    if filename.split('.')[-2].endswith('4'):
        grid=wave3 #dummy, to be consistent with the function output format.
        #exclude O2 absorption when measuring SNR 
        non_tel=np.where(wave>=telluric_region)[0]
        fl=d[non_tel]
        signal = np.median(fl)
        noise = 1.4826/np.sqrt(2)*np.median(np.abs(fl[1:] - fl[:-1]))
        snr = signal / noise
    
    print 'S/N is ' + str(snr)
    if snr <3 and filename.split('.')[-2].endswith('4')==False :
        sn=False
        sn_low.append(starno)
    grid_min=min(grid, key=lambda x:abs(x-wave[0]))
    grid_max=min(grid, key=lambda x:abs(x-wave[-1]))
    grid=grid[np.where(grid==grid_min)[0][0]:np.where(grid==grid_max)[0][0]]
    flux=np.interp(grid,wave,d)
    return flux,grid,sn,snr,d,wave

#to find the subpixel, fitting the correlation peaks
fitfunc2=lambda p,x: p[0]*x**2+p[1]*x+p[2]
errfunc2=lambda p,x,y: fitfunc2(p,x)-y

def find_rv (model_ctm,flux,wave,filename,ccd):
    
    '''
    inputs:
    1. model
    2. data, ctm normalised
    3. common wavelength gird
    4. full path to fits
    
    outputs:
    1. rv 
    2. max(correlation coeff) for this particular model 
    '''
    flux = flux-np.mean(flux) #get rid of the top hat like thingy
    h=pyfits.getheader(filename)
    if ccd==1:
        cdelt1=wave1[1]-wave1[0]
    if  ccd==2:
        cdelt1=wave2[1]-wave2[0]
    if ccd==3:
        cdelt1=wave3[1]-wave3[0]
#    if ccd==2 and filename=='/Users/jlin/guess/data/150426/combined/1504260006010011.fits':
#        pdb.set_trace()
    flux = flux/(np.sqrt(np.sum(flux**2))) #normalise the flux
    model_ctm1 = model_ctm- np.mean(model_ctm) #get rid of the top hat like thingy
    model_ctm1 /= np.sqrt(np.sum(model_ctm1**2)) #normalise the flux

    #window to smooth out the ends so it dont break when cross correlating yo 
    slope=np.arange(len(wave)*0.10)# creating a ramp using 10% of the spectrum on either side 
    window=np.hstack([slope,np.ones(len(wave)-len(slope)*2)*slope[-1],-1*slope+slope[-1]])
    window=window/max(window)
    model_ctm1=window*model_ctm1
    flux=window*flux

    #performing the cross correlation
    coeff=np.correlate(model_ctm1,flux,mode='same')

    max_index=np.where(coeff==max(coeff))[0][0]
    #fitting the cross correlation peak
    x_vals=np.array([max_index-1,max_index,max_index+1])
    d=x_vals[-1]-len(coeff)
    if d>=0:
        x_vals=np.array([max_index-2,max_index-1,max_index])  
    try:
        y_vals=coeff[x_vals]
    except IndexError:
        x_vals=np.array([max_index,max_index+1,max_index+2])
    y_vals=coeff[x_vals]
    p0,p1,p2=leastsq(errfunc2,[0,0,0],args=(x_vals,y_vals))[0]
    x=np.arange(min(x_vals),max(x_vals),0.1) 
    max_x=x[np.where(p0*x**2+p1*x+p2==max(p0*x**2+p1*x+p2))[0][0]]
    shift=max_x-len(model_ctm1)//2
    dlambda=cdelt1*shift
    velocity=-1*dlambda * 2.99792458e8 / (h['LAMBDAC']*1000)
    #print velocity,max(coeff)
    return velocity,max(coeff)

m1=pyfits.open(rvmodel1)
m2=pyfits.open(rvmodel2)
m3=pyfits.open(rvmodel3)

fitfunc3=lambda p,x: p[0]*x**4+p[1]*x**3+p[2]*x**2+p[3]*x+p[4]
errfunc3=lambda p,x,y: fitfunc3(p,x)-y

def split_list(alist, wanted_parts=1):
    length = len(alist)
    return [ alist[i*length // wanted_parts: (i+1)*length // wanted_parts] 
             for i in range(wanted_parts) ]

def find_rv2 (filename,ccd): #ccd=2 or 3 only!
    starno=filename.split('/')[-1][:-6]
    flux,grid,sn,snr,a,b=clean(filename[:-6]+str(ccd)+'.fits')
    if snr<3:
        print 'sn too low'
        return(0,0,0,0,0)

    if ccd==3:
        wavee=wave3
        mm=m3
        h_s=min(list(grid), key=lambda x:abs(x-6530))
        h_e=min(list(grid), key=lambda x:abs(x-6592))
        h_s_index=np.where(grid==h_s)[0][0]
        h_e_index=np.where(grid==h_e)[0][0]
        s = np.polynomial.Chebyshev.fit(np.hstack([grid[:h_s_index],grid[h_e_index:]]), np.hstack([flux[:h_s_index],flux[h_e_index:]]), deg=7)
        flux_ctm=flux/s(grid)

    if ccd==2:
        wavee=wave2
        mm=m2
        p0,p1,p2,p3,p4=leastsq(errfunc3,[0,0,0,0,0],args=(grid,flux))[0]
        fit=p0*grid**4+p1*grid**3+p2*grid**2+p3*grid+p4
        flux_ctm=flux/fit#(grid)
        
    flux_split=split_list(flux_ctm,4)
    flux_ctm1,flux_ctm2,flux_ctm3,flux_ctm4=flux_split[0],flux_split[1],flux_split[2],flux_split[3]

    grid_split=split_list(grid,4)
    grid1,grid2,grid3,grid4=grid_split[0],grid_split[1],grid_split[2],grid_split[3]

    
    model_s=np.where(wavee==grid[0])[0][0]
    model_e=np.where(wavee==grid[-1])[0][0]
    models=mm[0].data[:,model_s:model_e+1]

    model_s1=np.where(wavee==grid1[0])[0][0]
    model_e1=np.where(wavee==grid1[-1])[0][0]
    models1=mm[0].data[:,model_s1:model_e1+1]

    model_s2=np.where(wavee==grid2[0])[0][0]
    model_e2=np.where(wavee==grid2[-1])[0][0]
    models2=mm[0].data[:,model_s2:model_e2+1]

    model_s3=np.where(wavee==grid3[0])[0][0]
    model_e3=np.where(wavee==grid3[-1])[0][0]
    models3=mm[0].data[:,model_s3:model_e3+1]

    model_s4=np.where(wavee==grid4[0])[0][0]
    model_e4=np.where(wavee==grid4[-1])[0][0]
    models4=mm[0].data[:,model_s4:model_e4+1]    

    coeffs1,coeffs2,coeffs3,coeffs4,coeffs=[],[],[],[],[] 
    rvs1,rvs2,rvs3,rvs4,rvs=[],[],[],[],[]


    for j in range(len(models1))[0:-3]:
        rv1,coeff1=find_rv(models1[j],flux_ctm1,grid1,filename[:-6]+str(ccd)+'.fits',ccd)
        coeffs1.append(coeff1)
        rvs1.append(rv1)
        rv2,coeff2=find_rv(models2[j],flux_ctm2,grid2,filename[:-6]+str(ccd)+'.fits',ccd)
        coeffs2.append(coeff2)
        rvs2.append(rv2)
        rv3,coeff3=find_rv(models3[j],flux_ctm3,grid3,filename[:-6]+str(ccd)+'.fits',ccd)
        coeffs3.append(coeff3)
        rvs3.append(rv3)
        rv4,coeff4=find_rv(models4[j],flux_ctm4,grid4,filename[:-6]+str(ccd)+'.fits',ccd)
        coeffs4.append(coeff4)
        rvs4.append(rv4)

        rv,coeff=find_rv(models[j],flux_ctm,grid,filename[:-6]+str(ccd)+'.fits',ccd)
        coeffs.append(coeff)
        rvs.append(rv)


    thres=0.3
    good_coeff1=np.where(np.array(coeffs1)>thres)[0]
    good_rv1=np.array(rvs1)[good_coeff1]
    weights1=np.array(coeffs1)[good_coeff1]/float(sum(np.array(coeffs1)[good_coeff1]))

    good_coeff2=np.where(np.array(coeffs2)>thres)[0]
    good_rv2=np.array(rvs2)[good_coeff2]
    weights2=np.array(coeffs2)[good_coeff2]/float(sum(np.array(coeffs2)[good_coeff2]))

    good_coeff3=np.where(np.array(coeffs3)>thres)[0]
    good_rv3=np.array(rvs3)[good_coeff3]
    weights3=np.array(coeffs3)[good_coeff3]/float(sum(np.array(coeffs3)[good_coeff3]))

    good_coeff4=np.where(np.array(coeffs4)>thres)[0]
    good_rv4=np.array(rvs4)[good_coeff4]
    weights4=np.array(coeffs4)[good_coeff4]/float(sum(np.array(coeffs4)[good_coeff4]))

    good_coeff=np.where(np.array(coeffs)>thres)[0]
    good_rv=np.array(rvs)[good_coeff]
    weights=np.array(coeffs)[good_coeff]/float(sum(np.array(coeffs)[good_coeff]))

    combined_rv1=sum(weights1*good_rv1)
    combined_rv2=sum(weights2*good_rv2)
    combined_rv3=sum(weights3*good_rv3)
    combined_rv4=sum(weights4*good_rv4)
    combined_rv=sum(weights*good_rv)
    combined_rvs=[combined_rv1,combined_rv2,combined_rv3,combined_rv4]
    combined_rvs=[np.nan if x==0 else x for x in combined_rvs]
    
    return( combined_rvs[0], combined_rvs[1],combined_rvs[2],combined_rvs[3],combined_rv )

    
def run_rv (filename,folder):

    '''
    inputs:
    1. path to file, just one ccd, it finds the rest 2 automatically
    2. folder
    
    outputs:
    1. 3 rvs for 3 ccds 
    '''    
    starno=filename.split('/')[-1][:-6]
    rvss=[]
    snrs=[]
    #checks sn
    
    snr1=clean('%s/%s/combined/%s%s.fits'%(dir_str,folder,starno,'1'))[-3]#folder,date,ccd,galahic
    snr2=clean('%s/%s/combined/%s%s.fits'%(dir_str,folder,starno,'2'))[-3]#folder,date,ccd,galahic
    snr3=clean('%s/%s/combined/%s%s.fits'%(dir_str,folder,starno,'3'))[-3]#folder,date,ccd,galahic
    snr4=clean('%s/%s/combined/%s%s.fits'%(dir_str,folder,starno,'4'))[-3]

    snrs.append(snr1)
    snrs.append(snr2)
    snrs.append(snr3)
    snrs.append(snr4)

    if snr1<3:
        print 'yooo dwag sn are looowwwwww'
        return np.nan,snrs,1,np.nan

    for i in ['1','2','3']: #looping thru the ccds
        filename='%s/%s/combined/%s%s.fits'%(dir_str,folder,starno,i)
        flux,grid,sn,snr,a,b = clean(filename)                 
        if i =='1':
            wavee=wave1
            mm=m1
            #excluding H regions when fitting continua
            h_s=min(list(grid), key=lambda x:abs(x-4847))
            h_e=min(list(grid), key=lambda x:abs(x-4900))
            h_s_index=np.where(grid==h_s)[0][0]
            h_e_index=np.where(grid==h_e)[0][0]
            s = UnivariateSpline(np.hstack([grid[:h_s_index],grid[h_e_index:]]), np.hstack([flux[:h_s_index],flux[h_e_index:]]), s=1e15) #s=smoothing factor
            flux_ctm=flux/s(grid)

        if i == '2':
            wavee=wave2
            mm=m2
            p0,p1,p2,p3,p4=leastsq(errfunc3,[0,0,0,0,0],args=(grid,flux))[0]
            fit=p0*grid**4+p1*grid**3+p2*grid**2+p3*grid+p4
            flux_ctm=flux/fit#(grid)
          
        if i == '3':
            wavee=wave3
            mm=m3
            h_s=min(list(grid), key=lambda x:abs(x-6530))
            h_e=min(list(grid), key=lambda x:abs(x-6592))
            h_s_index=np.where(grid==h_s)[0][0]
            h_e_index=np.where(grid==h_e)[0][0]
            #s = UnivariateSpline(np.hstack([grid[:h_s_index],grid[h_e_index:]]), np.hstack([flux[:h_s_index],flux[h_e_index:]]), s=1e15)
            s = np.polynomial.Chebyshev.fit(np.hstack([grid[:h_s_index],grid[h_e_index:]]), np.hstack([flux[:h_s_index],flux[h_e_index:]]), deg=7)
            flux_ctm=flux/s(grid)

        model_s=np.where(wavee==grid[0])[0][0]
        model_e=np.where(wavee==grid[-1])[0][0]
        models=mm[0].data[:,model_s:model_e+1]

        coeffs=[] 
        rvs=[]

        #finding the rvs for all 15 models 
        for j in range(len(models))[0:-3]:
            rv,coeff=find_rv(models[j],flux_ctm,grid,filename,int(i))
            coeffs.append(coeff)
            rvs.append(rv)


        good_coeff=np.where(np.array(coeffs)>0.3)[0]
        #only take the rvs where the cross correlation is reasonable (ie cross correlation
        #coeff >0.3)
        good_rv=np.array(rvs)[good_coeff]

        weights=np.array(coeffs)[good_coeff]/float(sum(np.array(coeffs)[good_coeff]))
        #final rv for this ccd is the weighted sum of rvs from ~15 models, weighted by ccoeff
        print 'weighted rv ccd'+i, sum(weights*good_rv)
        if sum(weights*good_rv) ==0 and i=='2': #incase the cc is crap
            rvss.append(np.nan)
            continue
        if sum(weights*good_rv) ==0 and (i=='1' or i=='3'):
            print 'ohh bad weights'
            bad_weights.append(starno)
            return np.nan,snrs,np.nan,1
            break
        rvss.append(sum(weights*good_rv))
    
    if rvss != [] and len(rvss)==3 :        
        return rvss,snrs,0,0

def sorting1 (rvs): #rvs=[rv1,rv2,rv3]
    #this step excludes the outlier if it lies further than 2x the distance between the other two
    i=rvs
    i=[999999 if np.isnan(x)==True else x for x in i]

    good_stars=[] #[ rvs, sorted rvs, final rv]
    max_id,min_id,mid_id=i.index(max(i)),i.index(min(i)),i.index(heapq.nlargest(2,i)[1])
    if abs(i[max_id]-i[mid_id])/abs(i[mid_id]-i[min_id])>2:
            good_stars.append([i,[i[mid_id],i[min_id]],np.mean([i[mid_id],i[min_id]])])
            good_stars[0][0]=[np.nan if x==999999 else x for x in good_stars[0][0]]
            good_stars[0][1]=[np.nan if x==999999 else x for x in good_stars[0][1]]

            #print i,np.mean([i[mid_id],i[min_id]])
            return good_stars   
    if abs(i[mid_id]-i[min_id])/abs(i[max_id]-i[mid_id])>2:
            #print i, np.mean([i[mid_id],i[max_id]])
            good_stars.append([i,[i[mid_id],i[max_id]],np.mean([i[mid_id],i[max_id]])])
            good_stars[0][1]=[np.nan if x==999999 else x for x in good_stars[0][1]]
            good_stars[0][1]=[np.nan if x==999999 else x for x in good_stars[0][1]]
            return good_stars   
    good_stars.append([i,i,np.mean(i)])
    good_stars[0][0]=[np.nan if x==999999 else x for x in good_stars[0][0]]
    good_stars[0][1]=[np.nan if x==999999 else x for x in good_stars[0][1]]
    #print i,np.mean(i) 
    return good_stars     


rvs=[]
ids=[]
if guess_params[16] == 'all':
    files=glob.glob('%s/%s/combined/*1.fits'%(dir_str,folder))
else:
    files=glob.glob('%s/%s/combined/*1.fits'%(dir_str,folder))[:int(guess_params[16])]


f=open(saveloc,'w') #change store location
f.write('{0:<25} {1:<10} {2:<10} {3:<10} {4:<10} {5:<10} {6:<10} {7:<10} {8:<10} {9:<10} {10:<10} {11:<10} {12:<10} {13:<10}  {14:<10}  {15:<10} {16:<10} {17:<10} {18:<10} {19:<10} {20:<10} {21:<10}\n' .\
format('ccd1_filename','v_ccd1','v_ccd2','v_ccd3','v_final','v_sigma', 's/n_ccd1', 's/n_ccd2','s/n_ccd3', 's/n_ccd4','sn_low', 'bad_weights','v_ccd2_1','v_ccd2_2','v_ccd2_3','v_ccd2_4','v_ccd3_1','v_ccd3_2','v_ccd3_3','v_ccd3_4','std_v2','std_v3'))


for i in files:
    print i 
    rv=run_rv(i,folder)
    rv2=find_rv2(i,2)
    rv3=find_rv2(i,3)

    rv2_std=np.std(np.array(rv2[:-1])[np.where(np.isnan(np.array(rv2[:-1]))==False)])
    rv3_std=np.std(np.array(rv3[:-1])[np.where(np.isnan(np.array(rv3[:-1]))==False)])

    if rv[-1]==0:
        if np.isnan(rv[0][1])!=True and rv2[4]!=0:
            assert rv[0][1]==rv2[4]
            print 'ok'
        elif np.isnan(rv[0][2])!=True and rv3[4]!=0:
            assert rv[0][2]==rv3[4]
            print 'ok'
    star_id=i.split('/')[-1].split('.')[0]
    if rv[2]==1:#snr low
        f.write('{0:<25} {1:<10.3f} {2:<10.3f} {3:<10.3f} {4:<10.3f} {5:<10.3f} {6:<10.3f} {7:<10.3f} {8:<10.3f} {9:<10.3f} {10:<10} {11:<10} {12:<10.3f} {13:<10.3f}  {14:<10.3f}  {15:<10.3f} {16:<10.3f} {17:<10.3f} {18:<10.3f} {19:<10.3f} {20:<10.3f} {21:<10.3f}\n'.\
    format(star_id,np.nan,np.nan,np.nan,np.nan,np.nan,rv[1][0],rv[1][1],rv[1][2],rv[1][3],1,0,rv2[0],rv2[1],rv2[2],rv2[3],rv3[0],rv3[1],rv3[2],rv3[3],rv2_std,rv3_std))
        continue
    if rv[3]==1:
        f.write('{0:<25} {1:<10.3f} {2:<10.3f} {3:<10.3f} {4:<10.3f} {5:<10.3f} {6:<10.3f} {7:<10.3f} {8:<10.3f} {9:<10.3f} {10:<10} {11:<10} {12:<10.3f} {13:<10.3f}  {14:<10.3f}  {15:<10.3f} {16:<10.3f} {17:<10.3f} {18:<10.3f} {19:<10.3f} {20:<10.3f} {21:<10.3f}\n'.\
    format(star_id,np.nan,np.nan,np.nan,np.nan,np.nan,rv[1][0],rv[1][1],rv[1][2],rv[1][3],0,1,rv2[0],rv2[1],rv2[2],rv2[3],rv3[0],rv3[1],rv3[2],rv3[3],rv2_std,rv3_std))
        continue
      #bad_weights
    if np.isnan(rv[0][1])==True and rv[2]==0 and rv[3]==0:
        good_star=sorting1(rv[0])
        f.write('{0:<25} {1:<10.3f} {2:<10.3f} {3:<10.3f} {4:<10.3f} {5:<10.3f} {6:<10.3f} {7:<10.3f} {8:<10.3f} {9:<10.3f} {10:<10} {11:<10} {12:<10.3f} {13:<10.3f}  {14:<10.3f}  {15:<10.3f} {16:<10.3f} {17:<10.3f} {18:<10.3f} {19:<10.3f} {20:<10.3f} {21:<10.3f}\n'.\
    format(star_id,good_star[0][0][0],good_star[0][0][1],good_star[0][0][2],\
           good_star[0][2],np.std([good_star[0][0][0],good_star[0][0][-1]]),rv[1][0],rv[1][1],rv[1][2],rv[1][3],0,0,rv2[0],rv2[1],rv2[2],rv2[3],rv3[0],rv3[1],rv3[2],rv3[3],rv2_std,rv3_std)) ####change the std 
    else:
        good_star=sorting1(rv[0])
        f.write('{0:<25} {1:<10.3f} {2:<10.3f} {3:<10.3f} {4:<10.3f} {5:<10.3f} {6:<10.3f} {7:<10.3f} {8:<10.3f} {9:<10.3f} {10:<10} {11:<10} {12:<10.3f} {13:<10.3f}  {14:<10.3f}  {15:<10.3f} {16:<10.3f} {17:<10.3f} {18:<10.3f} {19:<10.3f} {20:<10.3f} {21:<10.3f}\n'.\
    format(star_id,good_star[0][0][0],good_star[0][0][1],good_star[0][0][2],\
           good_star[0][2],np.std(good_star[0][0]),rv[1][0],rv[1][1],rv[1][2],rv[1][3],0,0,rv2[0],rv2[1],rv2[2],rv2[3],rv3[0],rv3[1],rv3[2],rv3[3],rv2_std,rv3_std))  ####change the std 
    
f.close() 

