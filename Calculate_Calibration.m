function [coefficients, fval, exitflag, output] = Calculate_Calibration( ...
    phantom_img, spacing, phantom_radius, FOV_diameter, n_circle_steps, circle_width, ...
    options)
%% Calculate_Calibration - calculates the FDC
% Algorithm accepts phantom image, its spacing, radius and an optimset and
% returns the coefficients for FDC after calculating them.
% Note: Execution might take up to 15 min, depending on image size and
% options.
% 
% INPUTS:
%   * phantom_img : [k x n x m] array
%       OCT image data of spherical phantom. Dimensions are to be denoted
%       as [Z x X x Y], with Z being the dimension along the optical axis
%       (depth).
%   * spacing : [1 x 3] array
%       Spacing of the OCT imaging device ([Z x X x Y]).
%   * phantom_radius : float
%       REAL radius of spherical phantom.
%   * options : [1x1] struct
%       Optimization options for fminsearch-function. For more information,
%       search documentation on 'fminsearch' and 'optimset'.
%
% OUTPUTS:
%   * coefficients : struct
%       Coefficients for usage in Apply_Coefficients_To_Surface.
%   * fval : float
%       Final loss function value in the end of optimization process.
%   * exitflag : int
%       Exitflag for abortion criteria. For more information,
%       search documentation on 'fminsearch'.
%   * output : [1x1] struct
%       Information about the optimization process. For more information,
%       search documentation on 'fminsearch'.
%
% DEPENDENCIES / TOOLBOXES
%   * Apply_Coefficients_To_Surface.m
%   * The MATLAB Optimization Toolbox is required for the execution of
%   this file.
%
% Author: Maron Dolling
% Institute for Biomedical Optics - Universitaet zu Luebeck
% and
% Medical Laser Center Luebeck
% Email: m.dolling@uni-luebeck.de
% January 2023
%------------------------- START MAIN -------------------------------------

%% Initialize start values for field correction coefficients
% Note: coefficients will later be saved in a struct with notation in
% paper-correspondence
coefficients_initial = [
    0,1,0;      % Quadratic in x direction
    0,1,0;      % Quadratic in y direction
    1,0,0;      % z-dependent scaling of x- and y-direction
    0,0,0;      % |
    0,0,0;      % |
    0,0,0;      % | => Zernike surface distortion approximation
    0,0,0;      % |
    0,0,0;      % |
    ];


%% Prepare OCT data
phantom_img_surface = Surface_Detection_Phantom(phantom_img, spacing);


%% Start optimization process
[optim_out, fval, exitflag, output] = fminsearch( @(c) ...
        loss_function(phantom_img_surface, phantom_radius, c, FOV_diameter, ...
        n_circle_steps, circle_width), coefficients_initial, options);


%% Make coefficients a struct with corresponding notation
c.q10 = optim_out(1,1); c.q11 = optim_out(1,2); c.q12 = optim_out(1,3); 
c.q20 = optim_out(2,1); c.q21 = optim_out(2,2); c.q22 = optim_out(2,3);

c.s01 = optim_out(3,1); c.s02 = optim_out(3,2); c.s03 = optim_out(3,3);

c.c0 = optim_out(4,1); c.c1 = optim_out(4,2); c.c2 = optim_out(4,3);
c.c3 = optim_out(5,1); c.c4 = optim_out(5,2); c.c5 = optim_out(5,3);
c.c6 = optim_out(6,1); c.c7 = optim_out(6,2); c.c8 = optim_out(6,3);
c.c9 = optim_out(7,1); c.c10 = optim_out(7,2); c.c11 = optim_out(7,3);
c.c12 = optim_out(8,1); c.c13 = optim_out(8,2); c.c14 = optim_out(8,3);

coefficients = c;

end
%------------------------- END MAIN ---------------------------------------


%% Loss function for Optimization prcedure
function loss = loss_function(phantom_surface, phantom_radius, ...
        coefficients, FOV_diameter, n_circle_steps, circle_width)
    % Calculates the mean deviation of radii at certain angles around the
    % phantom from the expected.
    
    % Apply the coefficients to correct surface
    phantom_surface_corrected = Apply_Coefficients_To_Surface( ...
            phantom_surface, coefficients, FOV_diameter/2);
    
    x = phantom_surface_corrected(:,1);
    y = phantom_surface_corrected(:,2);
    z = phantom_surface_corrected(:,3);

    % Remove any nan-values from data
    ii = and(and(~isnan(x), ~isnan(y)), ~isnan(z));
    x = x(ii);
    y = y(ii);
    z = z(ii);

    % Find BFS and center the phantom in x-y-plane
    [center, ~] = ellipsoid_fit([x,y,z], 'xyz');
    x = x - center(1); 
    y = y - center(2); 
    z = z - center(3);
    
    % Get all the angles the circle is fittet at
    angles = linspace(-pi/2, pi/2 , n_circle_steps);
    % add small random, so angles are not the same in every iteration
    angles = angles + (rand * pi / n_circle_steps);
    % Handle angles that exceed the limit
    angles(angles > pi/2) = angles(angles > pi/2) - pi;
    angles = sort(angles);
    
    % Calculate the radius for every angle
    radii = nan(1, n_circle_steps);
    for i = 1:n_circle_steps
    
        % Calculate the two points l1, l2 that define the line that
        % encapsules all points with a distance smaller than "width"
        l1 = repmat([cos(angles(i)), sin(angles(i)), 0], length(x), 1);
        l2 = -l1;
        distances = point_to_line_distance([x,y,zeros(length(x),1)], l1, l2);
        
        % Get indices of points within cirlce widt
        ii = find(distances < circle_width / 2);
    
        % Skip sphere fit for very sparse surfaces
        if length(ii) < 4
            continue
        end
        
        % Calculate radius of circle at angle 
        [~, r] = ellipsoid_fit([x(ii), y(ii), z(ii)], 'xyz');
        radii(i) = r(1);

    end

    % Calculate loss as mean deviation of radii form expected
    loss = mean(abs(radii-phantom_radius));

end







