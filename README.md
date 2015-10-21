# GUESS
Mr. Bones' Wild Ride

Jane Lin 
u5027368@anu.edu.au
last updated 21/10/15

--------------------------------------CHANGE LOG---------------------------------------

- 21/10/15:
updated model wavelength gird from:

lscale1= np.loadtxt(plscale1)[:,0][200:-150]#4717.0471-4886.7621

lscale2= np.loadtxt(plscale2)[:,0][50:-150]#5651.2565- 5864.2118

lscale3= np.loadtxt(plscale3)[:,0][20:-150] #6483.1451-6729.374 

lscale4= np.loadtxt(plscale4)[:,0][50:-150]#7710.6577-7880.5470

to:
lscale1= np.loadtxt(plscale1)[:,0][354:3709]#4724.026-4876.022

lscale2= np.loadtxt(plscale2)[:,0][284:3723]#5564.05-5852.02

lscale3= np.loadtxt(plscale3)[:,0][337:3737] #6503.032-6716.26 

lscale4= np.loadtxt(plscale4)[:,0][50:-150]#7710.6577-7880.5470

(ambre4.5new.p)

- 17/10/15: 
added in warning for stars outside the model grid 
added new dir location in the guess_param.txt


--------------------------------------HOW TO RUN IT------------------------------------

- GUESS reads in the parameter file 'guess_param.txt', located in the same folder as guess_2dfdr.py 
guess_param.txt specifies where to keep the outputs and where the relevant files are
- GUESS accepts the following directory structure: 
?????/yymmdd/combined/?????.fits 
- it runs the whole /yymmdd/combined folder at once (or x number of stars specified in guess_param.txt)


--------------------------------------WHAT IT DOES--------------------------------------

- GUESS outputs RV, parameter estimates, RV corrected & normalised spectra (Normalised spectra 
have the suffix '_n' in it)


--------------------------------------KNOWN ISSUES--------------------------------------

- RV seems to fail for very hot stars


--------------------------------------NOTES---------------------------------------------

- only trust(ish) GUESS values for stars with 'out' flag=0 in the output
- normalised spectrum is not generated if the continuum/rv is not trustworthy 
- telluric regions (specified in guess_param.txt) are cropped off before RV correction
- relevant chapters in my thesis are 3 (RV) and 4 (parameters), prolonged reading may cause 
severe headache
