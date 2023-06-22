function [Htotal, tnew, td, ex, ey, ez, dDtmn, exm, eyn, B2x, B2y] = vpirOfRectangularArrayPistonlikeTransducersWScanner(a, b, c1, x, y, z, fs, N, M, kerf_x, kerf_y, delayLawEnabled, zf, xf, yf, F, D, discretization) % Main function
    a = abs(a); % Output absolute value of input
    b = abs(b); % Output absolute value of input
    z = abs(z);
    X = x;
    Y = y;
    Z = z;
	F = abs(F); % Output absolute value of input
    
% X-axis
m = 1:M; % The mth element in the array
sx = 2*a + kerf_x; % The distance between centroids (Pitch)
B2x = M*2*a + (M - 1)*kerf_x; % The total length of the array [m]

% Y-axis
n = 1:N; % The mth element in the array
sy = 2*b + kerf_y; % The distance between centroids (Pitch)
B2y = N*2*b + (N - 1)*kerf_y; % The total length of the array [m]

% Geometry parameters for the mth element of an array for use in multiple and single line
% source array beam models and for considering the far field response of the array. The distance
% em is measured from the origin (taken as the center of the array) to the centroid of the mth
% element while xn is the distance measured from this element centroid to the center of the nth
% segment of this element.

exm = zeros(length(n), length(m)); % Define the size of matrix
eyn = zeros(length(n), length(m)); % Define the size of matrix

for m = 1:M % Column(s)
    for n = 1:N % Row(s)
        exm(n, m) = sx*(2*m -1 - M)/2; % The x-location of the centroid of the mth element
        eyn(n, m) = sy*(2*n -1 - N)/2; % The y-location of the centroid of the mth element
    end
end

% Delay law
m=1:M;
n=1:N;
Mb = (M - 1)/2;
Nb = (N - 1)/2;
exm_bar = ((m - 1) - Mb)*sx;
eyn_bar = ((n - 1) - Nb)*sy;
Ux = X(:)./F;
Uy = Y(:)./F;
Uz = Z(:)./F;
rmn = sqrt((xf - exm).^2 + (yf - eyn).^2 + (zf).^2);

% NAN is not zero. Using isnan to check the NAN value and convert it to zero.

% https://mathworld.wolfram.com/Indeterminate.html
% https://www.mathworks.com/matlabcentral/answers/115782-multiplication-of-infinity-by-zero-in-matlab-calculation
% https://en.wikipedia.org/wiki/IEEE_754
% https://stackoverflow.com/questions/36016266/how-to-get-a-floating-point-infinity-that-when-multiplied-by-zero-gives-zero

% IEEE conventions state
% https://www.mathworks.com/matlabcentral/answers/4426-how-to-handle-nan-values

rmn(isnan(rmn))= 0; 
ux = (x - exm)./rmn;
uy = (y - eyn)./rmn;
uz = z./rmn;
ex = rmn.*ux; ex(isnan(ex))= 0;
ey = rmn.*uy; ey(isnan(ey))= 0;
ez = rmn.*uz; ez(isnan(ez))= 0;

switch(F)
    case (Inf)
        % Dtmn = (exm*sin(theta)*cos(phii) + eyn*sin(theta)*sin(phii))/c1;
        Dtmn = (exm*Ux + eyn*Uy)/c1;
        dDtmn = abs(min(min(Dtmn))) + Dtmn;
    otherwise        
% Considering steering and focusing of the array wave field to a point, P,
% in the surrounding medium specified by the local distance, F, and the
% steering angles (theta, phii). In this case, the distance, rmn, from an
% element to P is given from the geometry as:
        % rmn = sqrt((R*sin(theta)*cos(phii) - exm).^2 + (R*sin(theta)*sin(phii) - eyn).^2 + R^2*(cos(theta))^2);
        
        % Delay law    
        % Beam Steering Through a Planar Interface
        dDtmn = max(max(rmn/c1)) - rmn/c1;
end 

% Velocity potential impulse response of rectangular pistonlike transducers
h = cell(length(ey(:, 1)), length(ex(1, :)), length(ez(1, 1)));
t = cell(length(ey(:, 1)), length(ex(1, :)), length(ez(1, 1)));

for zz = 1:length(ez(1, 1))
    for yy = 1:length(ey(:, 1))
        for xx = 1:length(ex(1, :))
            [t_temp, h_temp] = summation(a, b, c1, ex(1, xx), ey(yy, 1), ez(zz, 1), fs, D, discretization);                           
            h{yy, xx, zz} = h_temp;
            t{yy, xx, zz} = t_temp;
        end
    end
end

% This loop creates the cell array for the time delays for an array with MxN elements.
td = cell(length(ey(:, 1)), length(ex(1, :)), length(ez(1, 1))); %// Initialize matrix
for zz = 1:length(ez(1, 1))
    for yy = 1:length(ey(:, 1))
        for xx = 1:length(ex(1, :))                    
            td{yy, xx, zz} = t{yy, xx, zz} + delayLawEnabled * dDtmn(yy, xx);
        end
    end
end 
%The sum total of velocity potential impulse response 
tmin = min(min(cellfun(@min,td)));
tmax = max(max(cellfun(@max,td)));
STEP = t{1, 1}(2) - t{1, 1}(1);
tnew= tmin:STEP:tmax;
fidx = zeros(size(td, 1), size(td, 2));  %// Initialize matrix
lidx = zeros(size(td, 1), size(td, 2));  %// Initialize matrix
for yy = 1:size(td, 1)
    for xx = 1:size(td, 2)
		% First position
		[~, fq] = min(abs(tnew - td{yy, xx}(1)));
		fidx(yy, xx) = fq;
		% Last position
		[~, lq] = min(abs(tnew - td{yy, xx}(end)));
		lidx(yy, xx) = lq;
    end
end
 
if max(max(cellfun(@length, td))) > length(tnew)
    r = abs(max(max(cellfun(@length, td))) - length(tnew));
    tnew = [tnew NaN*ones(1, r)];
end

Htemp = cell(size(td, 1), size(td, 2)); %// Initialize matrix
Htotal = zeros(length(tnew), 1);  %// Initialize matrix
for yy = 1:size(td, 1)
    for xx = 1:size(td, 2)
        % s = length(tnew) - abs(length(tnew) - fidx(yy, xx) - length(h{yy, xx}') - abs(length(tnew) - lidx(yy, xx)));
        s = length(tnew) - fidx(yy, xx) - length(h{yy, xx}') + 1;
        if s >= 1
            Htemp{yy, xx}  = [zeros(abs(fidx(yy, xx) - 1), 1)' h{yy, xx}' zeros(s, 1)'];
        elseif s == 0
            Htemp{yy, xx}  = [zeros(abs(fidx(yy, xx) - 1), 1)' h{yy, xx}'];
        elseif s < 0
            s = abs(length(tnew) - fidx(yy, xx) - length(h{yy, xx}'));
            Htemp{yy, xx}  = [zeros(abs(fidx(yy, xx) - s), 1)' h{yy, xx}'];
        else            
            if length(h{yy, xx}) == length(tnew)
                Htemp{yy, xx}  = h{yy, xx}';
            else
                Htemp{yy, xx}  = [zeros(abs(fidx(yy, xx) - 1), 1)' h{yy, xx}' zeros(lidx(yy, xx), 1)'];
            end 
        end
        Htotal = Htotal + Htemp{yy, xx}';        
    end
end
end