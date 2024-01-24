% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Velocity Potential Impulse Response Of Circular Pistonlike Transducers
%
% Name: vpirOfCircularPistonlikeTransducers.m
%
% Function: Based on the article "Ultrasonic Beam Structures in Fluid
% Media" from J.P. Weight on 1984, it calculates the pressure impulse
% response of a circular aperture assuming an ideal piston behaviour.
%
% SOFTWARE: MATLAB R2014b
%
% Author: Alexandre Rabelo
%
% Change Log: 
% Date:				Name:					Description:
% 2016/08/24 		Alexandre Rabelo		1.0 - Initial Version
% 2016/08/26 		Flávio Buiochi          1.0 - Reviewed, Tested and Approved
% 2019/05/05        Alexandre Rabelo        1.1 - Changed the variable N to fs
%
% To invoke this function, firstly run the script:
% 'loadingDataOfCircularPiston.m'
%
% Vectorized mathematical code is used to calculate.
%
% %%%%%% Required input
% x -- > Abscissa (X-axis) [m].
% y -- > Ordinate (Y-axis) [m].
% D -- > Diameter of US-transducer [m].
% c -- > The velocity of sound in the propagating medium [m/s].
% f -- > The operating frequency of US-transducer [Hz].
% fs -- > The sampling frequency or sampling rate, fs, is the average number of samples obtained in one second (samples per second), thus fs = 1/T;
%
% The velocity potential impulse response at a point Q in the field of an idealized piston source undergoing uniform motion and radiating into a lossless fluid medium.
%
function [t, h] = vpirOfCircularPistonlikeTransducersWScanner(x, y, D, c, fs) % Main function
    
    x = abs(x); % Output absolute value of input
    y = abs(y); % Output absolute value of input
    D = abs(D); % Output absolute value of input
    
    [r, R, t0, t1, t2, t] = calculating(x, y, D, c, fs);
    
    % The big omega is the angle of equidistant arc included on the surface.
    OMEGA = angleOfEquidistantArc(r, x, y, R, t0, t1, t2, t);
    
    % Velocity Potential Impulse Response 
    h = c*OMEGA/(2*pi);
        
    % Normalizing data to 0-1 range.
    % h = real(h/c);    
end

function [r, R, t0, t1, t2, t] = calculating(x, y, D, c, fs)
    t0 = x/c; % [s]
    R = D/2; % [m]
    t1 = c^-1 * sqrt((R-y)^2 + x^2); % [s]
    t2 = c^-1 * sqrt((R+y)^2 + x^2); % [s]
    
    % Frequency sample
    % Nyquist rate -- The sampling frequency should be at least twice the
    % highest frequency contained in the signal (fs >= 2*fc).
    
    delta_t = 1/fs;
 
	Nn = min([t0 t1 t2]);
	Nx = max([t0 t1 t2]);

    if (Nn - delta_t < 0)
        t = Nn : delta_t : Nx + delta_t; % [s]
    else
        t = Nn - delta_t : delta_t : Nx + delta_t; % [s]
    end

	r = c*t; % [m]  
end

% The big omega is the angle of equidistant arc included on the surface.
function OMEGA = angleOfEquidistantArc(r, x, y, R, t0, t1, t2, t)
    if (y < R)
        arc = InsideGeometrical(r, x, y, R, t0, t1, t2, t);
    elseif (y == R)
        arc = OnEdge(r, x, R, t0, t1, t2, t);
    else
        arc = OutsideGeometrical(r, x, y, R, t1, t2, t);
    end
    OMEGA = arc;
end

% Expressions for the angle of equidistant arc on the surface of a circular source.

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Region: Inside geometrical (Beam: y < R)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function arc = InsideGeometrical(r, x, y, R, t0, t1, t2, t)
    arc = zeros(size(t));
    arc(t < t0) = 0;
    arc((t0 <= t) & (t <= t1)) = 2*pi;
    arc((t1 < t) & (t <= t2)) = 2*acos((r((t1 < t) & (t <= t2)).^2 - x^2 + y^2 - R^2)./(2*y*sqrt(r((t1 < t) & (t <= t2)).^2 - x^2)));
    arc(t > t2) = 0;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Region: On Edge (Beam: y = R)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function arc = OnEdge(r, x, R, t0, t1, t2, t)
    arc = zeros(size(t));
    arc(t < t0) = 0;
    arc((t0 == t) & (t == t1)) = pi;
    arc((t1 < t) & (t <= t2)) = 2*acos(sqrt(r((t1 < t) & (t <= t2)).^2 - x^2)/(2*R));
    arc(t > t2) = 0;
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Region: Outside geometrical (Beam: y > R)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function arc = OutsideGeometrical(r, x, y, R, t1, t2, t)
    arc = zeros(size(t));
    arc(t <= t1) = 0;
    arc((t1 < t) & (t <= t2)) = 2*acos((r((t1 < t) & (t <= t2)).^2 - x^2 + y^2 - R^2)./(2*y*sqrt(r((t1 < t) & (t <= t2)).^2 - x^2)));
    arc(t > t2) = 0;
end

%
% Reference
%

% WEIGHT, J.P., "Ultrasonic Beam Structures in Fluid Media”, J. Acoust. Soc. Am., v.76, p.1184-1191, 1984.
% Support material for the vectorization instruction can be found at http://www.mathworks.com/help/matlab/matlab_prog/vectorization.html