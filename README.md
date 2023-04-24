# OCT_Distortion_Correction
Code for Publication of OCT Fan Distortion Correction

## Author
For questions, suggestions, etc. please don't hesitate to contact

Maron Dolling
Institute for Biomedical Optics - Universitaet zu Luebeck
and
Medical Laser Center Luebeck
Email: m.dolling@uni-luebeck.de


## Setup
The software was written and debugged in MATLAB R2022b. Additional requirements are
- the MATLAB Optimization tool box for the calculation of the correction
- the MATLAB image processing toolbox for the phantom surface detection.

If you do not have access to these, you might be able to adjut the algorithm and subsitute the Downhill-Simplex algorithm (fminsearch()) in "Calculate_Calibration.m" with another search strategy. Moreover, the surface detection algorithm may be adjusted to work without the bwconncomp() function. Adaptation of the phantom surface detection algorithm might be necessary anyway, depending on how the phantom image aqcuisition quality is with your spcific phantom and OCT setup.


## Scripts

### EXAMPLE.m
This script makes a full walkthrough through the calibration and coefficient application process. All necessary steps are shown here.
The software prerequisites that you are able to have your OCT data in a Z x X x Y - shaped MATLAB array. In the code, dummy data is inserted. WITHOUT SUBSITUTING THE DUMMY DATA, THE CODE WILL NOT WORK PROPERLY. OCT manufacturers like Thorlabs include MATLAB software that can export .oct-files to MATLAB arrays.

### Calculate_Calibration.m
Use this script to calculate the calibration. It returns a set of coefficients as a MATLAB struct, that has a notation accoding to the paper and can be applied using Apply_Coefficients_Surface.m.

### Apply_Coefficients_Surface.m
This is used to actually apply calibration coefficients onto surface data. It is also required to calculate the calibrtion.
Please make sure to only give extracted surface data complying with the scripts conventions (see MATLAB code documentation).

### Surface_Detection_Phantom.m
The surface detection is called within Calculate_Calibration.m. It may happen that your phantom or OCT setup delivers images that don't work well with this surface detection algorithm. In this case, you might have to adjust it.

### ellipsoid_fit.m
Fast sphere fit algorithm. Source: Yury Petrov - https://de.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit.

