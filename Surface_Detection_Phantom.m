function [SURFACE] = Surface_Detection_Phantom(IMG, spacing)
%% Surface_Detection_Phantom - detect surface of spherical phantom
% Algorithm for detection of the surface of a spherical phantom. Accepts
% the OCT data and voxel spacing and returns surface in cartesian
% coordinates.
% 
% INPUTS:
%   * IMG : [k x n x m] array
%       OCT image data of spherical phantom. Dimensions are to be denoted
%       as [Z x X x Y], with Z being the dimension along the optical axis
%       (depth).
%   * spacing : [1 x 3] array
%       Spacing of the OCT imaging device ([Z x X x Y]).
%
% OUTPUTS:
%   * SURFACE : [n*m x 3] array
%       Surface data with cartesian x-, y- and z-cooridnates, represented 
%       with the matrix's columns.
%
% DEPENDENCIES / TOOLBOXES
%   - MATLAB image processing toolbox
%
% Author: Maron Dolling
% Institute for Biomedical Optics - Universitaet zu Luebeck
% and
% Medical Laser Center Luebeck
% Email: m.dolling@uni-luebeck.de
% January 2023
%--------------------------------------------------------------------------

if ndims(IMG) ~= 3
    error("The phantom OCT image should be provided as an array of " + ...
        "size Z x X x Y.");
end
if length(spacing) ~= 3
    error("Please provide spacing data that corresponds to the image " + ...
        "dimensions");
end

% Read sizes
Nx = size(IMG,2);
Ny = size(IMG,3);

% Read spacing
dz = spacing(1);
dx = spacing(2);
dy = spacing(3);

% Get indices of max-values in every a-scan.
[~, ii] = max(IMG, [], 1);
ii = squeeze(ii);        
ii = round(filloutliers(ii, NaN));

% Create a logical image of size of input an assign surface points as
% "true" values
surface_logical = false(size(IMG));
for x = 1:Nx
    for y = 1:Ny
        if ~isnan(ii(x,y))
            surface_logical(ii(x,y), x, y) = true;
        end
    end
end 

% Detect biggest connected component to avoid speckle noise
new = false(size(surface_logical));
CC = bwconncomp(surface_logical);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~, idx] = max(numPixels);
new(CC.PixelIdxList{idx}) = 1;
surface_logical = new;

% Get indices again, and remove outliers of BCC due to effects close to
% edges of image or reflection etc.
[~, ii] = max(surface_logical, [], 1);
ii = squeeze(ii);        
ii(ii < 10) = NaN;

% Make a meshgrid to transform surface to cartesian coordiantes, and center
% the surface data around (0,0)
[x, y] = meshgrid(1:Nx, 1:Ny);
z = ii(:) * dz;
x = (x(:) - Nx/2) * dx;
y = (y(:) - Ny/2) * dy;

SURFACE = [x, y, z];

end 
