-------------------------------------------
GALAH ULTRA ENIGMATIC SPECTROSCOPIC SCRIPT

Jane Lin 
u5027368@anu.edu.au
-------------------------------------------

A glorified random number generator for getting RV, parameter estimates, RV corrected & normalised spectra for 2dfdr reduced GALAH data

Quick Start Guide:
- Run getmodels.py in models/ to download the model grid (alternatively, you can get it here https://www.dropbox.com/s/zy4d7k0fmwbbi4e/ambre4.5new.p)
- Edit the parameter file 'guess_param.txt' to specify where you want to keep everything
- Run guess_2dfdr.py
- ???????
- PROFIT! 


Notes:
- GUESS accepts the following directory structure: 
?????/yymmdd/combined/?????.fits 
- It runs the whole /yymmdd/combined folder at once (or x number of stars specified in guess_param.txt)
- Only trust(ish) GUESS values for stars with 'out' flag=0 in the output
- Normalised spectrum ('????_n') is not generated if the continuum/rv is not trustworthy 
- Telluric regions (specified in guess_param.txt) are cropped off before RV correction
- GUESS relevant chapters in my thesis are 3 (RV) and 4 (parameters) (prolonged reading may cause severe headache)


Known Issues:
- RV seems to fail for very hot stars
