%% EXAMPLE - Example of FDC calculation
% 
% This script shows an example on how to use the FDC algorithm. It is
% assumed that the user can provide the OCT data as a MATLAB array already.
% Thorlabs provides MATLAB software for this purpose.
% 
%
% Author: Maron Dolling
% Institute for Biomedical Optics - Universitaet zu Luebeck
% and
% Medical Laser Center Luebeck
% Email: m.dolling@uni-luebeck.de
% January 2023
%--------------------------------------------------------------------------

%% Input data
phantom_img = ones([512,512,512]); % this is your spherical phantom image
spacing = [0.00615, 0.01563, 0.01563]; % voxel spacing of OCT (z * x * y)
phantom_radius = 8.5; % real radius of phantom used
n_circle_steps = 30; % Number of steps the semicricle should be fitted at
circle_width = 0.1; % width of semicircle to be fitted

%--------------------------------------------------------------------------
%% Use phantom to calculate FDC
% set options for optimization process
% for more details, see 'optimset' documentation
options = optimset('Display', 'iter', 'TolFun',...
    1e-5, 'TolX', 1e-05, 'PlotFcns', @optimplotfval, 'MaxIter', 500);

% execute optimization
% for more details, see 'fminsearch' documentation
% Hint: the surface is extracted in Calculate_Calibration()
[coefficients, fval, exitflag, output] = Calculate_Calibration( ...
    phantom_img, spacing, phantom_radius, n_circle_steps, circle_width, ...
    options);

% 'coefficients' now carries the parameters that can be used in
% 'Apply_Coefficients_To_Surface' to correct any extracted surface (e.g.
% cornea) by changing its (x,y,z)-coordinates.

%% Using the phantom surface as an example
SURFACE = Surface_Detection_Phantom(phantom_img, spacing);
SURFACE_corrected = Apply_Coefficients_To_Surface(SURFACE, coefficients);
