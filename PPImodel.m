function losv = PPImodel(b,y,r)
% =========================================================================
% AUTHOR: Nicola Bodini - Sept 2016
% =========================================================================
% PPImodel returns an array of line-of-sight velocity (losv) measurements
% given a set of wake model parameters (b), an array of lateral distances
% (y), and a range gate (r).
%
% m = length(b) specifies which model is used:
%       m = 2: no-wake model
%       m = 14: 4 wakes model
% =========================================================================
% INPUT      DESCRIPTION                           UNITS     SIZE
% -------------------------------------------------------------------------
% b          wake model parameters                           [14x1]
%            b(1) = turbine yaw angle              [rad]
%            b(2) = ambient wind speed             [m/s]
%            b(3:6) = velocity deficit amplitudes  [m/s]
%            b(7:10) = wakes peak locations        [D]
%            b(10:14) = wakes widths               [D]
%
% y          lateral distances                     [D]       [nx1]
%
% r          range gate                            [D]       [1x1]
% =========================================================================
% OUTPUT     DESCRIPTION                           UNITS     SIZE
% -------------------------------------------------------------------------
% losv       Scanning Lidar losv measurements      [m/s]     [nx1]
% =========================================================================

% number of model parameters
m = length(b);

% turbine yaw angle
phi = b(1);
% ambient wind speed
u = b(2);

if m > 2 % 4 wakes model
    % velocity deficit amplitudes
    a(1) = b(3);
    a(2) = b(4);
    a(3) = b(5);
    a(4) = b(6);
    % wakes peak locations
    mu(1) = b(7);
    mu(2) = b(8);
    mu(3) = b(9);
    mu(4) = b(10);
    % wakes widths / 4
    s(1) = b(11);
    s(2) = b(12);
    s(3) = b(13);
    s(4) = b(14);
end


% actual flowfield
switch m
    case 2 % no wake case
        uActual = u;
    case 14 % 4 wakes -> Gaussian model
        uActual = u - a(1)*exp(-(y - mu(1)).^2/(2*s(1)^2))...
                    - a(2)*exp(-(y - mu(2)).^2/(2*s(2)^2))...  
                    - a(3)*exp(-(y - mu(3)).^2/(2*s(3)^2))...
                    - a(4)*exp(-(y - mu(4)).^2/(2*s(4)^2));
end

% OUTPUT: modeled line-of-sight velocity
% losv = u * cos(theta-phi)
losv = (sqrt(r^2-y.^2)/r*cos(phi) + y/r*sin(phi)).*uActual;