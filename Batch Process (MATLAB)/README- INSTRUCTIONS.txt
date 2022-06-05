Thank you for downloading IVEN and the batch version of IVEN in partiuclar. 

We would like to encourage you to use the full version of IVEN complete with GUIs for cell classifcation and threshold method selection. IVEN should not be used as a black box method and we encourage you to review any results output from IVEN thoroughly and in biological context. 
However, lots of data analsysi can be time consuming and we appreciate it is sometimes useful for a 'quick and dirty' intial approach. 

This batch version of IVEN uses no GUIs and simply runs through the selected files with given parameters. 

Be sure to define these parameters in the run_batch_IVEN.m file. (Namely the method for thresholding and associated parameters, as well as the shrink factor used within the convex hull algorithm to classify inside versus outside cells). 

If you have any questions, please refer to the standard tutorials or get in touch!

Currently we have only developed the batch processing method in MATLAB as this method appeared to be the most intensive when analysing large numbers of cells. 

Post-analysis figure inspection NOTE: 
-If you want to inspect the .fig files after running the batch process version of IVEN, please do the following. 
openfig('filename.fig') %make sure you're in the correct directory of the file
f=gcf;
set(f,'visible','on')   %this step is required because when the fig files were created, their visibility was set to 'off' to prevent lost of windows. 
