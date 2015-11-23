"""
Galah
Ultra 
Enigmatic
Spectroscopic
Script

last updated: 21-10-15

@author: Jane Lin (u5027368@anu.edu.au)

*use at yo own risk

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
    
    if filename.split('.')[0].endswith('1'):
        grid=wave1
    if filename.split('.')[0].endswith('2'):
        grid=wave2
    if filename.split('.')[0].endswith('3'):
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

    if filename.split('.')[0].endswith('4'):
        grid=wave3 #dummy, to be consistent with the function output format.
        #exclude O2 absorption when measuring SNR 
        non_tel=np.where(wave>=7680)[0]
        fl=d[non_tel]
        signal = np.median(fl)
        noise = 1.4826/np.sqrt(2)*np.median(np.abs(fl[1:] - fl[:-1]))
        snr = signal / noise
    
    print 'S/N is ' + str(snr)
    if snr <3:
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

def find_rv (model_ctm,flux,wave,filename):
    
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

    flux -= np.mean(flux) #get rid of the top hat like thingy
    h=pyfits.getheader(filename)
    if filename.split('.')[0].endswith('1'):
        cdelt1=wave1[1]-wave1[0]
    if  filename.split('.')[0].endswith('2'):
        cdelt1=wave2[1]-wave2[0]
    if filename.split('.')[0].endswith('3'):
        cdelt1=wave3[1]-wave3[0]
 
    flux /= np.sqrt(np.sum(flux**2)) #normalise the flux
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

def run_rv (filename,folder):

    '''
    inputs:
    1. path to file, just one ccd, it finds the rest 2 automatically
    2. folder
    
    outputs:
    1. 3 rvs for 3 ccds 
    '''    

    #starno=filename.split('_')[-1].split('.')[0]#galahic
    starno=filename.split('/')[-1][:-6]
    #date=filename.split('/')[-1].split('_')[0][:-1]
    #date=filename[:-1] #not really date
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
        return 999,snrs,1,999

    for i in ['1','2','3']: #looping thru the ccds
        filename='%s/%s/combined/%s%s.fits'%(dir_str,folder,starno,i)
        flux,grid,sn,snr,a,b = clean(filename)        
        snrs.append(snr)
            
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
            flux_ctm=flux/fit            
        if i == '3':
            wavee=wave3
            mm=m3
            h_s=min(list(grid), key=lambda x:abs(x-6530))
            h_e=min(list(grid), key=lambda x:abs(x-6592))
            h_s_index=np.where(grid==h_s)[0][0]
            h_e_index=np.where(grid==h_e)[0][0]
            s = UnivariateSpline(np.hstack([grid[:h_s_index],grid[h_e_index:]]), np.hstack([flux[:h_s_index],flux[h_e_index:]]), s=1e15)
            flux_ctm=flux/s(grid)
        model_s=np.where(wavee==grid[0])[0][0]
        model_e=np.where(wavee==grid[-1])[0][0]
        models=mm[0].data[:,model_s:model_e+1]

        coeffs=[] 
        rvs=[]
        #finding the rvs for all 15 models 
        for j in range(len(models))[0:-3]:
            rv,coeff=find_rv(models[j],flux_ctm,grid,filename)
            coeffs.append(coeff)
            rvs.append(rv)
        good_coeff=np.where(np.array(coeffs)>0.3)[0]
        #only take the rvs where the cross correlation is reasonable (ie cross correlation
        #coeff >0.3)
        good_rv=np.array(rvs)[good_coeff]
        #pdb.set_trace()
        weights=np.array(coeffs)[good_coeff]/float(sum(np.array(coeffs)[good_coeff]))
        #final rv for this ccd is the weighted sum of rvs from ~15 models, weighted by ccoeff
        print 'weighted rv', sum(weights*good_rv)
        if sum(weights*good_rv) ==0 and i=='2': #incase the cc is crap
            rvss.append(999)
            continue
        if sum(weights*good_rv) ==0 and (i=='1' or i=='3'):
            print 'ohh bad weights'
            bad_weights.append(starno)
            return 999,snrs,999,1
            break
        rvss.append(sum(weights*good_rv))
    
    if rvss != [] and len(rvss)==3 :        
        return rvss,snrs,0,0

def sorting1 (rvs): #rvs=[rv1,rv2,rv3]
    #this step excludes the outlier if it lies further than 2x the distance between the other two
    i=rvs
    good_stars=[] #[ rvs, sorted rvs, final rv]
    max_id,min_id,mid_id=i.index(max(i)),i.index(min(i)),i.index(heapq.nlargest(2,i)[1])
    if abs(i[max_id]-i[mid_id])/abs(i[mid_id]-i[min_id])>2:
            good_stars.append([i,[i[mid_id],i[min_id]],np.mean([i[mid_id],i[min_id]])])
            #print i,np.mean([i[mid_id],i[min_id]])
            return good_stars   
    if abs(i[mid_id]-i[min_id])/abs(i[max_id]-i[mid_id])>2:
            #print i, np.mean([i[mid_id],i[max_id]])
            good_stars.append([i,[i[mid_id],i[max_id]],np.mean([i[mid_id],i[max_id]])])
            return good_stars   
    good_stars.append([i,i,np.mean(i)])
    #print i,np.mean(i) 
    return good_stars     


rvs=[]
ids=[]
if guess_params[16] == 'all':
    files=glob.glob('%s/%s/combined/*1.fits'%(dir_str,folder))
else:
    files=glob.glob('%s/%s/combined/*1.fits'%(dir_str,folder))[:int(guess_params[16])]


#donotwant=glob.glob('/priv/miner3/galah/testsample_05282015/%s/%s*'%(folder,folder.split('/')[-1]))
f=open(saveloc,'w') #change store location
f.write('{0:<25} {1:<10} {2:<10} {3:<10} {4:<10} {5:<10} {6:<10} {7:<10} {8:<10} {9:<10} {10:<10} {11:<10}\n' .\
format('ccd1_filename','v_ccd1','v_ccd2','v_ccd3','v_final','v_sigma', 's/n_ccd1', 's/n_ccd2','s/n_ccd3', 's/n_ccd4','sn_low', 'bad_weights'))


for i in files:
    print i 
#    if i in donotwant:
#	continue
    rv=run_rv(i,folder)
    
    star_id=i.split('/')[-1].split('.')[0]
    if rv[2]==1:#snr low
        f.write('{0:<25} {1:<10.1f} {2:<10.1f} {3:<10.1f} {4:<10.1f} {5:<10.1f} {6:<10.1f} {7:<10.1f} {8:<10.1f} {9:<10.1f} {10:<10} {11:<10}\n'.\
    format(star_id,999,999,999,999,999,rv[1][0],rv[1][1],rv[1][2],rv[1][3],1,0))
        continue
    if rv[3]==1:
        f.write('{0:<25} {1:<10.1f} {2:<10.1f} {3:<10.1f} {4:<10.1f} {5:<10.1f} {6:<10.1f} {7:<10.1f} {8:<10.1f} {9:<10.1f} {10:<10} {11:<10}\n'.\
    format(star_id,999,999,999,999,999,rv[1][0],rv[1][1],rv[1][2],rv[1][3],0,1))
        continue
      #bad_weights
    if rv[0][1]==999 and rv[2]==0 and rv[3]==0:
        good_star=sorting1(rv[0])
        f.write('{0:<25} {1:<10.1f} {2:<10.1f} {3:<10.1f} {4:<10.1f} {5:<10.1f} {6:<10.1f} {7:<10.1f} {8:<10.1f} {9:<10.1f} {10:<10} {11:<10}\n'.\
    format(star_id,good_star[0][0][0],good_star[0][0][1],good_star[0][0][2],\
    good_star[0][2],np.std([good_star[0][0][0],good_star[0][0][-1]]),rv[1][0],rv[1][1],rv[1][2],rv[1][3],0,0))
    else:
        good_star=sorting1(rv[0])
        f.write('{0:<25} {1:<10.1f} {2:<10.1f} {3:<10.1f} {4:<10.1f} {5:<10.1f} {6:<10.1f} {7:<10.1f} {8:<10.1f} {9:<10.1f} {10:<10} {11:<10}\n'.\
    format(star_id,good_star[0][0][0],good_star[0][0][1],good_star[0][0][2],\
    good_star[0][2],np.std(good_star[0][0]),rv[1][0],rv[1][1],rv[1][2],rv[1][3],0,0))
    
f.close() 

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return idx

#For generating the rv shifted spectrum, median filtered and uninterpolated
rvs=np.genfromtxt(saveloc, skiprows=1,usecols=[0,4],dtype='S')
for i in rvs: #writes all the rv corrected spectra
    #date=i[0].split('_')[0][:-1]
    #galahic=i[0].split('_')[1]
    date=i[0][:-1]
    if int(float(i[-1]))==999:
        continue
    rv=float(i[1])*-1000
    print 'shifting ', i 
    for j in ['1','2','3']:
        filename='%s/%s/combined/%s%s.fits'%(dir_str,folder,date,j)
        cleann=clean(filename)
        flux,wave=cleann[-2],cleann[-1]
        rest_wave=wave*(1+rv/2.99792458e8)
        with open('%s/%s%s.txt' %(folderloc,date,j), 'wb') as f:
            writer = csv.writer(f)
            writer.writerows(izip(rest_wave,flux))
    
    # For shifting the 4th CCD 
    filename='%s/%s/combined/%s4.fits'%(dir_str,folder,date)
    h = pyfits.getheader(filename)
    x = np.arange(h['NAXIS1']) + 1 - h['CRPIX1']
    wave = h['CDELT1']*x + h['CRVAL1']    
    hdu=pyfits.open(filename)
    d=hdu[0].data
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
    signal = np.median(d)
    noise = 1.4826/np.sqrt(2)*np.median(np.abs(d[1:] - d[:-1]))
    med_d = sig.medfilt(d,5) #5 pixels window
    w = np.where(abs(med_d - d) > sigma * noise)[0]
    d[w]=med_d[w]    
    telluric=find_nearest(wave,telluric_region)
    wave,d=wave[telluric:],d[telluric:]
    rest_wave=wave*(1+rv/2.99792458e8)
    with open('%s%s4.txt' %(folderloc,date), 'wb') as f:
        writer = csv.writer(f)
        writer.writerows(izip(rest_wave,d))



rv_input=saveloc
rv_corrected=folderloc

lscale1= np.loadtxt(plscale1)[:,0][354:3709]#4724.026-4876.022
lscale2= np.loadtxt(plscale2)[:,0][284:3723]#5564.05-5852.02
lscale3= np.loadtxt(plscale3)[:,0][337:3737] #6503.032-6716.26 
lscale4= np.loadtxt(plscale4)[:,0][50:-150]#7710.6577-7880.5470

lscale_all=np.hstack([lscale1,lscale2,lscale3])

stars=np.genfromtxt(rv_input, skiprows=1, dtype='S')
star_no=stars[:,0]
v1=np.array([float(i) for i in stars[:,1]])
v2=np.array([float(i) for i in stars[:,2]])
v3=np.array([float(i) for i in stars[:,3]])
vf=np.array([float(i) for i in stars[:,4]])
v_sig=np.array([float(i) for i in stars[:,5]])
sn1=np.array([float(i) for i in stars[:,6]])
sn2=np.array([float(i) for i in stars[:,7]])
sn3=np.array([float(i) for i in stars[:,8]])
sn4=np.array([float(i) for i in stars[:,9]])
sn_low=np.array([int(i) for i in stars[:,10]])
bad_weights=np.array([int(i) for i in stars[:,11]])

hermes_ctm=np.loadtxt(continuum_reg)

#uncomment this for output file
f=open(param_output,'w')
f.write('{0:<18} {1:<10} {2:<10} {3:<10} {4:<10} {5:<10} {6:<10} {7:<10} {8:<10} {9:<10} {10:<10} {11:<10} {12:<10} {13:<10} {14:<10} {15:<10} {16:<10}\n' .\
format('ccd1_filename','v_ccd1','v_ccd2','v_ccd3','v_final','v_sigma','s/n_ccd1', 's/n_ccd2','s/n_ccd3','s/n_ccd4','sn_low','bad_weights','ctm','teff','logg','feh','out'))


datas=[]#contains ctm normalized data
ctm=[]
for i in star_no: 
   index=np.where(star_no==i)[0][0]
   if int(vf[index])== 999:
       ctm.append(1)
       print 'bad rv, skipped'
       continue
   data=[]
   #date,galahic=i.split('_')[0][:-1],i.split('_')[1]
   date=i[:-1]
   print i
   fluxintt=[]
   for j in ['1','2','3','4']:
       bad=False
       fname=rv_corrected+date+j+'.txt'
       if fname.split('.')[-2].endswith('1'):
           lscale= lscale1
           sf=1e6
       if fname.split('.')[-2].endswith('2'):
           lscale= lscale2
           sf=1e6
       if fname.split('.')[-2].endswith('3'):
           lscale= lscale3
           sf=1e8
       if fname.split('.')[-2].endswith('4'):
           lscale= lscale4
           sf=1e9
       data1=np.loadtxt(fname,delimiter=',')
       flux=data1[:,1]
       raw_wave=data1[:,0]
       if j!='4' and (data1[:,0][0]>lscale[0] or data1[:,0][-1]<lscale[-1]):
       #if len(data1)<len(lscale):
           bad=True
           print 'ccd%s, star wl within the model wl, skipped'%j
           break
       fit_y=[]
       fit_x=[]
       for ss in range(len(raw_wave)):
           for kk in hermes_ctm:
             if raw_wave[ss]< kk[1] and raw_wave[ss]> kk[0]:
                 fit_y.append(flux[ss])
                 fit_x.append(raw_wave[ss])
       s = UnivariateSpline(fit_x, fit_y, s=sf)
       flux_raw_norm=flux/s(raw_wave)
       flux_int=np.interp(lscale,raw_wave,flux)
       flux_norm=flux_int/s(lscale)

       if j=='2' and len(np.where(flux_norm>1.6)[0])>5:
                   s = UnivariateSpline(fit_x, fit_y, s=1e7)
                   flux_int=np.interp(lscale,raw_wave,flux)
                   flux_norm=flux_int/s(lscale)
                   flux_raw_norm=flux/s(raw_wave)


       if np.isnan(flux_norm)[0] and  np.isnan(flux_norm)[-1] and \
np.isnan(flux_norm)[1000]:
           print 'fitting gone wrong!'
           bad=True
           break
       print len(flux_norm)
       if j!= '4':
       	   data.append(flux_norm)#each ccd is normalized individually
       	   if bad:
           	print 'woah, fast star! or something went wrong! skipped'

       with open(rv_corrected+date+j+'_n.txt', 'wb') as nom:
           writer = csv.writer(nom)
           writer.writerows(izip(raw_wave,flux_raw_norm))
   if bad:
       ctm.append(1)
       continue
   assert len(np.hstack([data[0],data[1],data[2]]))==10194
   datas.append(np.hstack([data[0],data[1],data[2]]))
   ctm.append(0)

datas=np.array(datas)
model=pickle.load(open(param_grid,'rb'))
models=model[0]
mparams=model[1]

teff_m=np.array([s[0] for s in mparams])
logg_m=np.array([s[1] for s in mparams])
feh_m=np.array([s[2] for s in mparams])

print 'neighbor search'

def distance_matrix(data,model):
    data=np.transpose(data)
    distance = np.tile(np.sum(data**2,axis=0),model.shape[0]).reshape(model.shape[0],data.shape[1]) - \
        2*np.dot(model,data) + \
        np.sum(model**2,axis=1).repeat(data.shape[1]).reshape(model.shape[0],data.shape[1])
    return distance #each col is for one star

distance=distance_matrix(datas,models)

nk=10 #number of neighbors here, can change
teff_ordered=sorted(list(set(teff_m)))
feh_ordered=sorted(list(set(feh_m)))

#%%
def neighbor(distance_matrix,data,i):#i=range(len(datas)) 
    neighs=[]
    ind=[]
    for j in range(nk):	
        index=np.where(distance_matrix[:,i]==np.array(sorted(distance_matrix[:,i])[j]))
        neighs.append(models[index][0])
        ind.append(index[0][0])
    if len(set(teff_m[ind]))==1:
        if teff_m[ind[0]]==2500:
            dteff=min(distance_matrix[:,i][np.where(teff_m==2700)[0]])
        if teff_m[ind[0]]==8000:
            dteff=min(distance_matrix[:,i][np.where(teff_m==7750)[0]])
        else:
            teffp,teffm=teff_ordered[np.where( teff_ordered==teff_m[ind[0]])[0][0]+1],\
teff_ordered[np.where( teff_ordered==teff_m[ind[0]])[0][0]-1]
            dteffp=min(distance_matrix[:,i][np.where(teff_m==teffp)[0]])
            dteffm=min(distance_matrix[:,i][np.where(teff_m==teffm)[0]])
            if dteffp>dteffm:
                dteff=dteffm
            else:
                dteff=dteffp
        neighs[-1]=models[np.where(distance_matrix[:,i]==dteff)[0][0]]
        ind[-1]=np.where(distance_matrix[:,i]==dteff)[0][0]
    if len(set(logg_m[ind]))==1:
        if logg_m[ind[0]]==-0.5:
            dlogg=min(distance_matrix[:,i][np.where(logg_m==0)[0]])
        if logg_m[ind[0]]==5.5:
            dlogg=min(distance_matrix[:,i][np.where(logg_m==5)[0]])
        else:
            loggp,loggm=logg_m[ind[0]]+0.5,logg_m[ind[0]]-0.5
            #pdb.set_trace()
            dloggp=min(distance_matrix[:,i][np.where(logg_m==loggp)[0]])
            dloggm=min(distance_matrix[:,i][np.where(logg_m==loggm)[0]])
            if dloggp>dloggm:
                dlogg=dloggm
            else:
                dlogg=dloggp
        neighs.append(models[np.where(distance_matrix[:,i]==dlogg)[0][0]])
        ind.append(np.where(distance_matrix[:,i]==dlogg)[0][0])

    if len(set(feh_m[ind]))==1:
        if feh_m[ind[0]]==-5:
            dfeh=min(distance_matrix[:,i][np.where(feh_m==-4)[0]])
        if feh_m[ind[0]]==1:
            dfeh=min(distance_matrix[:,i][np.where(feh_m==0.75)[0]])
        else:
            fehp,fehm=feh_ordered[np.where( feh_ordered==feh_m[ind[0]])[0][0]+1],\
feh_ordered[np.where( feh_ordered==feh_m[ind[0]])[0][0]-1]
            dfehp=min(distance_matrix[:,i][np.where(feh_m==fehp)[0]])
            dfehm=min(distance_matrix[:,i][np.where(feh_m==fehm)[0]])
            if dfehp>dfehm:
                dfeh=dfehm
            else:
                dfeh=dfehp
        neighs.append(models[np.where(distance_matrix[:,i]==dfeh)[0][0]])
        ind.append(np.where(distance_matrix[:,i]==dfeh)[0][0])

    return ind,neighs

print 'neighbor search over'

new_points=datas


def decompo(data,inds,neigh,i):#i=range(len(datas))
    print 'decomp ',i
    coeff=(np.array(neigh)).T #dp[0][i] is [[neighb1 coeff],[neighb2 coeff]...], want to take the transpose of this
    decomp=lstsq(coeff,data)
    nil=np.zeros(len(models))
    for kk in range(len(inds)):
        nil[inds[kk]]=decomp[0][kk]
    return nil

temps=[]
logg=[]
feh=[]
ww=[]
trolol=[]
for i in range(len(datas)):
    inds,neighs=neighbor(distance,datas[i],i)
    nil=(decompo(datas[i],inds,neighs,i))
    trolol.append(inds)
    w=nil/sum(nil)
    temps.append(np.dot(w.T,teff_m))
    logg.append(np.dot(w.T,logg_m))
    feh.append(np.dot(w.T,feh_m))
    ww.append(w)
ww=np.array(ww)

temps=np.array(temps)
logg=np.array(logg)
feh=np.array(feh)


ctm=np.array(ctm)
param_array=np.zeros([len(ctm),4])

param_array[:,0][np.where(ctm==0)[0]]=temps
param_array[:,1][np.where(ctm==0)[0]]=logg
param_array[:,2][np.where(ctm==0)[0]]=feh
param_array[:,0][np.where(ctm==1)[0]]=9999
param_array[:,1][np.where(ctm==1)[0]]=9999
param_array[:,2][np.where(ctm==1)[0]]=9999

#final flag, only use the params with out=0
dd=np.zeros((len(param_array),1))
dd[np.where(param_array[:,0]>8000)]=1
dd[np.where(param_array[:,0]<2500)]=1
dd[np.where(param_array[:,1]<-0.5)]=1
dd[np.where(param_array[:,1]>5.5)]=1
dd[np.where(param_array[:,2]<-5)]=1
dd[np.where(param_array[:,2]>1)]=1

#uncomment this for output file
for i in range(len(ctm)):
   f.write('{0:<18} {1:<10.1f} {2:<10.1f} {3:<10.1f} {4:<10.1f} {5:<10.1f} {6:<10.1f} {7:<10.1f} {8:<10.1f} {9:<10} {10:<10} {11:<10} {12:<10.1f} {13:<10.3f} {14:<10.3f} {15:<10} {16:<10}\n'.\
   format(star_no[i],v1[i],v2[i],v3[i],vf[i],v_sig[i],sn1[i],sn2[i],sn3[i],sn4[i],sn_low[i],\
   bad_weights[i],ctm[i],param_array[:,0][i],param_array[:,1][i],param_array[:,2][i],dd[i][0]))
f.close()
