        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Velocity Potential Impulse Response Of Rectangular Pistonlike Transducers
%
% Name: vpirOfRectangularPistonlikeTransducers.m
%
% Function: Based on the article "Diffraction impulse response of
% rectangular transducers" from Jose Luis San Emeterio and Luis G. Ullate 
% on 1992. In this aforementioned paper, the classical time-domain impulse
% response approach is used. The complexity introduced by the geometrical 
% discontinuities of the rectangular aperture is analysed in detail. A new
% compacting methodology is proposed and applied to the derivation of a
% closed-form expression for the diffraction impulse response of a
% uniformly vibrating rectangular transducer. Analytical expressions
% of h(x, t) are presented, dividing the positive quadrant into four
% geometrical regions:
% 
% Region I: x>=a and y>=b;
% Region II: x<a and y>=b;
% Region III: x>=a and y<b; 
% Region IV: x<a and y<b. 
%
% This new solution that was purposed by the authors provides the value
% of h(x, t) directly in the time domain for different boundary
% conditions, without requiring superposition methods.
%
% SOFTWARE: MATLAB R2014b
%
% Author: Alexandre Rabelo
%
% Change Log: 
% Date:				Name:					Description:
% 2016/08/24 		Alexandre Rabelo		1.0 - Initial Version.
% 2016/08/26 		Flávio Buiochi          1.0 - Reviewed, Tested and Approved.
% 2017/09/07 		Alexandre Rabelo        1.1 - Changed the size of the vector t.
% 2017/10/08        Alexandre Rabelo        1.2 - Changed the function inverseCircular.
% 2018/04/24        Marcelo Matuda          1.3 - Restored the function inverseCircular to the original one and added min(..., 1) for alpha(i)
% 2019/05/05        Alexandre Rabelo        1.4 - Changed the variable N to fs

%
% To invoke this function, firstly run the script: 
% 'loadingDataOfRetactangularPiston.m'
%
% Vectorized mathematical code is used to calculate.
%
% %%%%%% Required input
% a -- > Half length in [m].
% b -- > Half width/height in [m].
% c -- > The velocity of sound in the propagating medium [m/s].
% x -- > Abscissa (X-axis) [m].
% y -- > Ordinate (Y-axis) [m].
% z -- > Applicate (Z-axis) [m].
% f -- > The operating frequency of US-transducer [Hz]
% fs -- > The sampling frequency or sampling rate, fs, is the average number of samples obtained in one second (samples per second), thus fs = 1/T;
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% The velocity potential impulse response at a point P in the field of an idealized piston source undergoing uniform motion and radiating into a lossless fluid medium.
%
function [t, h] = vpirOfRectangularPistonlikeTransducers(a, b, c, x, y, z, fs) % Main function

    a = abs(a); % Output absolute value of input
    b = abs(b); % Output absolute value of input
    x = abs(x); % Output absolute value of input
    y = abs(y); % Output absolute value of input
    z = abs(z); % Output absolute value of input
    
    [sigma, d1, d2, d3, d4, TauA, TauB, TauC, TauD, Tau0, TauS1, TauS2, TauS3, TauS4, Taum, TauM, t] = calculating(a, b, c, x, y, z, fs);
	[alpha1, alpha2, alpha3, alpha4, alpha1_, alpha2_, alpha3_, alpha4_] = inverseCircular(sigma, d1, d2, d3, d4, Tau0, TauS1, TauS2, TauS3, TauS4, t);    
    
    % The big omega is the angle of equidistant arc included on the surface.
    OMEGA = geometricalRegion(alpha1, alpha2, alpha3, alpha4, alpha1_, alpha2_, alpha3_, alpha4_, TauA, TauB, TauC, TauD, Tau0, TauS1, TauS2, Taum, TauM, t, a, b, x, y);
    
    % Velocity Potential Impulse Response 
    h = c*OMEGA/(2*pi);
	
    % Normalizing data to 0-1 range.
    % h = h/c;
end

function [sigma, d1, d2, d3, d4, TauA, TauB, TauC, TauD, Tau0, TauS1, TauS2, TauS3, TauS4, Taum, TauM, t] = calculating(a, b, c, x, y, z, fs)
    
    % Frequency sample
    % Nyquist rate -- The sampling frequency should be at least twice the
    % highest frequency contained in the signal (fs >= 2*fc).
    
    STEP = 1/fs;
	
    % The sides of the rectangle
	d1 = x - a;
	d3 = x + a;
	d2 = y - b;
	d4 = y + b;

	TauA = sqrt(d1^2 + d2^2 + z^2)/c; % [s]
	TauB = sqrt(d2^2 + d3^2 + z^2)/c; % [s]
	TauC = sqrt(d1^2 + d4^2 + z^2)/c; % [s]
	TauD = sqrt(d3^2 + d4^2 + z^2)/c; % [s]
	Tau0 = z/c; % [s]
    
    % The times \tauSi are the active arcs that are tangent to sides Si (edge discontinuities)
	TauS1 = sqrt(d1^2 + z^2)/c; % [s]
	TauS2 = sqrt(d2^2 + z^2)/c; % [s]
	TauS3 = sqrt(d3^2 + z^2)/c; % [s]
	TauS4 = sqrt(d4^2 + z^2)/c; % [s]
	
	Taum = min(TauB, TauC);
	TauM = max(TauB, TauC);

	Nn = min([TauA TauB TauC TauD Tau0 TauS1 TauS2 TauS3 TauS4]);
	Nx = max([TauA TauB TauC TauD Tau0 TauS1 TauS2 TauS3 TauS4]);
    t = Nn - 10*STEP: STEP : Nx + 10*STEP; % [s]
  	r = c*t; % [m]
    
    % Radius
	sigma = sqrt(r.^2 - z^2); % r.^2 - z^2 must be greater than or equal to zero
end

%
% From an analytical and a computational point of views, each \alphai is the principal value of the inverse circular functions in the inverse sine [asin(')].
%

%function [alpha1, alpha2, alpha3, alpha4, alpha1_, alpha2_, alpha3_, alpha4_] = inverseCircular(sigma, d1, d2, d3, d4, Tau0, TauS1, TauS2, TauS3, TauS4, t)
%	alpha1(t >= TauS1) = asin(min(d1./sigma(t >= TauS1), 1));
%	alpha2(t >= TauS2) = asin(min(d2./sigma(t >= TauS2), 1));
%	alpha3(t >= TauS3) = asin(min(d3./sigma(t >= TauS3), 1));
%	alpha4(t >= TauS4) = asin(min(d4./sigma(t >= TauS4), 1));
%
%    % Alpha bar
%	alpha1_(t >= Tau0) = sign(d1)*asin(min(abs(d1)./sigma(t >= Tau0), 1));
%	alpha2_(t >= Tau0) = sign(d2)*asin(min(abs(d2)./sigma(t >= Tau0), 1));
%	alpha3_(t >= Tau0) = sign(d3)*asin(min(abs(d3)./sigma(t >= Tau0), 1));
%	alpha4_(t >= Tau0) = sign(d4)*asin(min(abs(d4)./sigma(t >= Tau0), 1));
%
%	alpha1_(t < TauS1) = sign(d1)*pi/2;
%	alpha2_(t < TauS2) = sign(d2)*pi/2;
%	alpha3_(t < TauS3) = sign(d3)*pi/2;
%	alpha4_(t < TauS4) = sign(d4)*pi/2;
%
%	alpha1_(t >= TauS1) = alpha1(t >= TauS1);
%	alpha2_(t >= TauS2) = alpha2(t >= TauS2);
%	alpha3_(t >= TauS3) = alpha3(t >= TauS3);
%	alpha4_(t >= TauS4) = alpha4(t >= TauS4);
%end

function [alpha1, alpha2, alpha3, alpha4, alpha1_, alpha2_, alpha3_, alpha4_] = inverseCircular(sigma, d1, d2, d3, d4, Tau0, TauS1, TauS2, TauS3, TauS4, t)
	alpha1(t >= TauS1) = asin(min(d1./sigma(t >= TauS1), 1));
	alpha2(t >= TauS2) = asin(min(d2./sigma(t >= TauS2), 1));
	alpha3(t >= TauS3) = asin(min(d3./sigma(t >= TauS3), 1));
	alpha4(t >= TauS4) = asin(min(d4./sigma(t >= TauS4), 1));

     % Alpha bar
     alpha1_(t > Tau0) = sign(d1)*asin(min(abs(d1)./sigma(t > Tau0), 1));
     alpha2_(t > Tau0) = sign(d2)*asin(min(abs(d2)./sigma(t > Tau0), 1));
     alpha3_(t > Tau0) = sign(d3)*asin(min(abs(d3)./sigma(t > Tau0), 1));
     alpha4_(t > Tau0) = sign(d4)*asin(min(abs(d4)./sigma(t > Tau0), 1));

     alpha1_(t < TauS1) = sign(d1)*pi/2;
     alpha2_(t < TauS2) = sign(d2)*pi/2;
     alpha3_(t < TauS3) = sign(d3)*pi/2;
     alpha4_(t < TauS4) = sign(d4)*pi/2;

      alpha1_(t > TauS1) = alpha1(t > TauS1);
      alpha2_(t > TauS2) = alpha2(t > TauS2);
      alpha3_(t > TauS3) = alpha3(t > TauS3);
      alpha4_(t > TauS4) = alpha4(t > TauS4);
end

%
% The expressions given in this function were extracted from Table I which is based on the rigid baffle boundary condition.
%

function OMEGA = geometricalRegion(alpha1, alpha2, alpha3, alpha4, alpha1_, alpha2_, alpha3_, alpha4_, TauA, TauB, TauC, TauD, Tau0, TauS1, TauS2, Taum, TauM, t, a, b, x, y)
    if (x >= a) && (y >= b) % Region I
        Taumin = TauA;
        arc(t <= Taumin) = 0;
        arc((Taumin < t) & (t < TauA)) = 0;
        arc((TauA <= t) & (t < Taum)) = pi/2 - alpha1((TauA <= t) & (t < Taum)) - alpha2((TauA <= t) & (t < Taum));
        if TauB <= TauC
            arc((Taum <= t) & (t < TauM)) = -alpha1((Taum <= t) & (t < TauM)) + alpha3((Taum <= t) & (t < TauM));
        else %TauB > TauC
            arc((Taum <= t) & (t < TauM)) = -alpha2((Taum <= t) & (t < TauM)) + alpha4((Taum <= t) & (t < TauM));
        end
        arc((TauM <= t) & (t <= TauD)) = -pi/2 + alpha3((TauM <= t) & (t <= TauD)) + alpha4((TauM <= t) & (t <= TauD));
        arc(t > TauD) = 0;
    elseif ((x < a) && (y >= b)) % Region II
        Taumin = TauS2;
        arc(t <= Taumin) = 0;
        arc((Taumin < t) & (t < TauA)) = pi - 2*alpha2((Taumin < t) & (t < TauA));
        arc((TauA <= t) & (t < Taum)) = pi/2 - alpha1((TauA <= t) & (t < Taum)) - alpha2((TauA <= t) & (t < Taum));
        if TauB <= TauC
            arc((Taum <= t) & (t < TauM)) = -pi - alpha1((Taum <= t) & (t < TauM)) + alpha3((Taum <= t) & (t < TauM)) + 2*alpha4_((Taum <= t) & (t < TauM));
        else  %TauB > TauC
            arc((Taum <= t) & (t < TauM)) = 0;
        end
        arc((TauM <= t) & (t <= TauD)) = -pi/2 + alpha3((TauM <= t) & (t <= TauD)) + alpha4((TauM <= t) & (t <= TauD));
        arc(t > TauD) = 0;
    elseif ((x >= a) && (y < b)) % Region III
        Taumin = TauS1;
        arc(t <= Taumin) = 0;
        arc((Taumin < t) & (t < TauA)) = 2*alpha3_((Taumin < t) & (t < TauA)) - 2*alpha1((Taumin < t) & (t < TauA));
        arc((TauA <= t) & (t < Taum)) = -pi/2 - alpha1((TauA <= t) & (t < Taum)) - alpha2((TauA <= t) & (t < Taum)) + 2*alpha3_((TauA <= t) & (t < Taum));
        if TauB <= TauC
            arc((Taum <= t) & (t < TauM)) = -alpha1((Taum <= t) & (t < TauM)) + alpha3((Taum <= t) & (t < TauM));
        else %TauB > TauC
            arc((Taum <= t) & (t < TauM)) = -pi - alpha2((Taum <= t) & (t < TauM)) + 2*alpha3_((Taum <= t) & (t < TauM)) + alpha4((Taum <= t) & (t < TauM));
        end
        arc((TauM <= t) & (t <= TauD)) = -pi/2 + alpha3((TauM <= t) & (t <= TauD)) + alpha4((TauM <= t) & (t <= TauD));
        arc(t > TauD) = 0;
    else % Region IV -- (x < a) && (y < b)
        Taumin = Tau0;
        arc(t <= Taumin) = 0;
        arc((Taumin < t) & (t < TauA)) = -2*pi - 2*alpha1_((Taumin < t) & (t < TauA)) - 2*alpha2_((Taumin < t) & (t < TauA)) + 2*alpha3_((Taumin < t) & (t < TauA)) + 2*alpha4_((Taumin < t) & (t < TauA));
        arc((TauA <= t) & (t < Taum)) = -3*pi/2 - alpha1((TauA <= t) & (t < Taum)) - alpha2((TauA <= t) & (t < Taum)) + 2*alpha3_((TauA <= t) & (t < Taum)) + 2*alpha4_((TauA <= t) & (t < Taum));
        if TauB <= TauC
            arc((Taum <= t) & (t < TauM)) = -pi - alpha1((Taum <= t) & (t < TauM)) + alpha3((Taum <= t) & (t < TauM)) + 2*alpha4_((Taum <= t) & (t < TauM));
        else %TauB > TauC
            arc((Taum <= t) & (t < TauM)) = -pi - alpha2((Taum <= t) & (t < TauM)) + 2*alpha3_((Taum <= t) & (t < TauM)) + alpha4((Taum <= t) & (t < TauM));
        end
        arc((TauM <= t) & (t <= TauD)) = -pi/2 + alpha3((TauM <= t) & (t <= TauD)) + alpha4((TauM <= t) & (t <= TauD));
        arc(t > TauD) = 0;
    end
    OMEGA = arc;
end

%
% Reference
%

% Support material for the vectorization instruction can be found at http://www.mathworks.com/help/matlab/matlab_prog/vectorization.html