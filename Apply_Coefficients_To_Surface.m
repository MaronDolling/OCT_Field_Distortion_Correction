function Y = Apply_Coefficients_To_Surface(X, c)
%% Apply_Coefficients_To_Surface -  correct an extracted surface
% Accepts surface data and coefficients for field distortion correction
% (FDC) and returns corrected surface coordinates.
% 
% INPUTS:
%   * X : [n x 3] array
%       Surface data with cartesian x-, y- and z-cooridnates, represented 
%       with the matrix's columns.
%   * c : struct
%       Coefficients for correction. Notation must be in correspondence to
%       the one in publication. Should be generated using
%       Calculate_Caibration().
%
% OUTPUTS:
%   * Y : [n x 3] array
%       Corrected surface data with cartesian x-, y- and z-cooridnates, 
%       represented with the matrix's columns.
%
% Author: Maron Dolling
% Institute for Biomedical Optics - University of Lübeck
% and
% Medical Laser Center Luebeck
% Email: m.dolling@uni-luebeck.de
% January 2023
%--------------------------------------------------------------------------

if size(X,2) ~= 3
    error("Please provide surface as cartesian " + ...
        "coordinates of shape n x 3.");
end

if ~isa(c, "struct")
    if size(c,1) ~= 8 || size(c,2) ~= 3
        error(['Invalid coefficients. Please check if given coefficients are provided' ...
            'by ''Calculate_Calibration'', or have correct notation. ' ...
            'If error occurs during execution of ''Calculate_Calibration'',' ...
            'the start values for the coefficients probably were changed ' ...
            'to invalid values.'])
    end
    try
        c_in = c;
        clear c;
        
        c.q10 = c_in(1); c.q11 = c_in(2); c.q12 = c_in(3); c.q13 = c_in(4); 
        c.q20 = c_in(5); c.q21 = c_in(6); c.q22 = c_in(7); c.q23 = c_in(8);
        
        c.s01 = c_in(9); c.s02 = c_in(10); c.s03 = c_in(11);
        
        c.c0 = c_in(12); c.c1 = c_in(13); c.c2 = c_in(14);
        c.c3 = c_in(15); c.c4 = c_in(16); c.c5 = c_in(17);
        c.c6 = c_in(18); c.c7 = c_in(19); c.c8 = c_in(20);
        c.c9 = c_in(21); c.c10 = c_in(22); c.c11 = c_in(23);
        c.c12 = c_in(24); c.c13 = c_in(25); c.c14 = c_in(26);

    catch
        error(['The coefficients could not be applied onto the input ' ...
            'surface. Please check if given coefficients are provided' ...
            'by ''Calculate_Calibration'', or have correct notation. ' ...
            'If error occurs during execution of ''Calculate_Calibration'',' ...
            'the start values for the coefficients probably were changed ' ...
            'to invalid values.'])
    end
end

x = X(:,1);
y = X(:,2);
z = X(:,3);

x = (c.q10 + c.q11*x + c.q12*x.^2) .* (c.s01 + c.s02*z + c.s03*z.^2 + y * c.y1); % x'(x)
y = (c.q20 + c.q21*y + c.q22*y.^2).* (c.s01 + c.s02*z + c.s03*z.^2 + x * c.x1); % y'(x)

[th, r] = cart2pol(x,y);

z = z + ...
        c.c0 + ...                 % offset
        c.c1 * (2*r.*sin(th)) + ... % ytilt
        c.c2 * (2*r.*cos(th)) + ... % xtilt
        c.c3 * sqrt(6) * (r.^2.*sin(2*th)) + ... % asti obl
        c.c4 * sqrt(3) * (2*r.^2-1) + ... % defocus
        c.c5 * sqrt(6) * (r.^2.*cos(2*th)) + ... % asti vert
        c.c6 * sqrt(8) * (r.^3.*sin(3*th)) + ... % trefoil vert
        c.c7 * sqrt(8) * (3*r.^3-2*r).*sin(th) + ... % vertical coma
        c.c8 * sqrt(8) * (3*r.^3-2*r).*sin(th) + ... % horizontal coma
        c.c9 * sqrt(8) * (r.^3.*cos(3*th)) + ... % trefoil obl
        c.c10 * sqrt(10) * (r.^4.*sin(4*th)) + ... % quadrafoil obl
        c.c11 * sqrt(10) * (4*r.^4 - 3*r.^2) .* sin(2*th) + ... % obl sec asti
        c.c12 * sqrt(5) * (6*r.^4 - 6*r.^2 + 1) + ... % primary spherical
        c.c13 * sqrt(10) * (4*r.^4 - 3*r.^2) .* cos(2*th) + ... % vert sec asti
        c.c14 * sqrt(10) * (r.^4.*cos(4*th)); % quadrafoil vert

Y = [x(:), y(:), z(:)];

end
