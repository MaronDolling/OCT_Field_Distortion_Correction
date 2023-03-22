function [BFS_radius, R_steep, A_steep, R_flat, A_flat, MSE, exitflag, ... 
    output] = Surface_Curvature_Assessment(SURFACE)
%% Corneal_Curvatue_Assessment - Get BFS and astigmatism of SURFACE
% Algorithm accepts phantom image, its spacing, radius and an optimset and
% returns the coefficients for FDC after calculating them.
% Note: Execution might take up to 15 min, depending on image size and
% options.
% 
% INPUTS:
%   * SURFACE : [n x 3] array
%       Surface data with cartesian x-, y- and z-cooridnates, represented 
%       with the matrix's columns.
%
% OUTPUTS:
%   * BFS_radius : float
%       Radius of best fitting sphere in the surface.
%   * R_steep : float
%       Steep astigmatism radius.
%   * A_steep : float
%       Steep astigmatism angle.
%   * R_flat : float
%       Flat astigmatism radius.
%   * A_flat : float
%       Flat astigmatism angle.
%   * MSE : float
%       Resdiual mean square error of surface fit.
%   * exitflag : int
%       Code for fit abortion cirterion. For more information, see
%       documentation on 'fminsearch'.
%   * output : struct
%       Different fit output metrics. For more information, see
%       documentation on 'fminsearch'.
%   
%
% DEPENDENCIES / TOOLBOXES
%   - MATLAB Optimization Toolbox
%   - point_to_line_distance.m
%   - ellipsoid_fit.m
%
% Author: Maron Dolling
% Institute for Biomedical Optics - Universitaet zu Luebeck
% and
% Medical Laser Center Luebeck
% Email: m.dolling@uni-luebeck.de
% March 2023
%--------------------------------------------------------------------------

% First: Remove NaN Entries from surface
ii = isnan(sum(SURFACE, 2));
SURFACE(ii, :) = [];

%% Calculate BFS
[~, R, ~] = ellipsoid_fit(SURFACE, 'xyz');
BFS_radius = R(1);

%% Calculate Zernike Fit
% Set options for the Zernike Fit
options = optimset('TolFun', 1e-10, 'TolX', 1e-10, ...
    'MaxIter', 1000, 'Display', 'none');

% set coefficients for optimization
coeffs_start = rand(1,6);

% set norm radius via surface coordinates
norm_radius = ceil( ...
    sqrt( ...
        max(SURFACE(:,1).^2 + SURFACE(:,2).^2, [], 'all') ...
        ) ...
    );

% Minimize fit function
[coeffs_optim, MSE, exitflag, output] = fminsearch( @(coeffs) ...
        Zernike_fit_loss(coeffs, SURFACE, norm_radius), ...
        coeffs_start, options);

% Get extremal angle
extremal_angle1 = 0.5 * atan(coeffs_optim(6) / coeffs_optim(4));
extremal_angle2 = wrapToPi(extremal_angle1 + pi/2);

% Get radii at extremal angles
R1 = get_radius_at_angle(SURFACE, extremal_angle1);
R2 = get_radius_at_angle(SURFACE, extremal_angle2);

% Check if steep or flat
if R1 < R2
    R_steep = R1;
    A_steep = extremal_angle1;
    R_flat = R2;
    A_flat = extremal_angle2;
else
    R_steep = R2;
    A_steep = extremal_angle2;
    R_flat = R1;
    A_flat = extremal_angle1;
end

end


%--------------------------------------------------------------------------


%% Loss Function for Zernike-Fit algorithm
function [loss] = Zernike_fit_loss(c, SURFACE, norm_radius)
    % gather x,y,z-coordinates
    x = SURFACE(:,1);
    y = SURFACE(:,2);
    z = SURFACE(:,3);

    % make polar and normalize
    [th, r] = cart2pol(x,y);
    r = r./norm_radius;

    % calc Zernike
    z = z + ...
        c(1) + ...                 % offset
        c(2) * (2*r.*sin(th)) + ... % ytilt
        c(3) * (2*r.*cos(th)) + ... % xtilt
        c(4) * sqrt(6) * (r.^2.*sin(2*th)) + ... % asti obl
        c(5) * sqrt(3) * (2*r.^2-1) + ... % defocus
        c(6) * sqrt(6) * (r.^2.*cos(2*th)); % asti vert

    % calculate deviation of fitted z to input z
    loss = immse( z, SURFACE(:,3) );
end

function [R] = get_radius_at_angle(SURFACE, alpha)
    
    width = 0.1; % standard width value

    % read surface data
    x = SURFACE(:,1); y = SURFACE(:,2); z = SURFACE(:,3);

    % Span line at angle and calculate distances of points to it
    l1 = repmat([cos(alpha), sin(alpha), 0], length(x), 1);
    l2 = -l1;
    distances = point_to_line_distance([x,y,zeros(length(x),1)], l1, l2);

    % get indices for valid points and calculate the semicircle radius
    ii = find(distances < width / 2);

    % Fit sphere into semicircle 
    [~, R, ~] = ellipsoid_fit([x(ii),y(ii),z(ii)], 'xyz');
    R = R(1);

end

    