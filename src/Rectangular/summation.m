a = 8e-3/2;           % a -- > Half length in [m].
b = 19e-3/2;          % b -- > Half width/height in [m].
c1 = 1500;            % c1 -- > The velocity of sound in the propagating medium [m/s].
fs = 16e6;            % Sampling frequency [Hz]
xc = 0e-3;            % xc --> x-coordinate of the center of the hydrophone
yc= 0e-3;             % yc --> y-coordinate of the center of the hydrophone
zc = 40e-3;           % zc --> z-coordinate for translation to the coordinate xy-plane
D = 6e-4;             % D --> Diameter of the hydrophone
discretization = 5;   % discretization --> Value to divide the geometry into finite elements to prepare for analysis

[xscanned, yscanned, zscanned, Nscanned] = scanner(xc, yc, zc, D, discretization);

% Velocity potential impulse response of rectangular pistonlike transducers
h = cell(1, Nscanned); %// Initialize matrix
t = cell(1, Nscanned); %// Initialize matrix

for index = 1:Nscanned
  [t_temp, h_temp] = vpirOfRectangularPistonlikeTransducers(a, b, c1, xscanned(index, 1), yscanned(index, 1), zscanned(index, 1), fs);
  h{index} = h_temp;
  t{index} = t_temp;
end

%The sum total of velocity potential impulse response 
tmin = min(cellfun(@min, t));
tmax = max(cellfun(@max, t));
STEP = t{1,1}(2) - t{1,1}(1);
tnew= tmin:STEP:tmax;
fidx = zeros(size(t, 1), size(t, 2));  %// Initialize matrix
lidx = zeros(size(t, 1), size(t, 2));  %// Initialize matrix
for yy = 1:size(t, 1)
    for xx = 1:size(t, 2)
		% First position
		[~, fq] = min(abs(tnew - t{yy, xx}(1)));
		fidx(yy, xx) = fq;
		% Last position
		[~, lq] = min(abs(tnew - t{yy, xx}(end)));
		lidx(yy, xx) = lq;
    end
end
 
if max(cellfun(@length, t)) > length(tnew)
    r = abs(max(cellfun(@length, t)) - length(tnew));
    tnew = [tnew NaN*ones(1, r)];
end

Htemp = cell(size(t, 1), size(t, 2)); %// Initialize matrix
Htotal = zeros(length(tnew), 1);  %// Initialize matrix
for yy = 1:size(t, 1)
    for xx = 1:size(t, 2)
        if abs(length(tnew) - fidx(yy, xx) - length(h{yy, xx}') - abs(length(tnew) - lidx(yy, xx))) > 1
            s = abs(length(tnew) - fidx(yy, xx) - length(h{yy, xx}') - abs(length(tnew) - lidx(yy, xx)));
            if length(h{yy, xx}) == length(tnew)
                Htemp{yy, xx}  = h{yy, xx}';
            else
                Htemp{yy, xx}  = [zeros(fidx(yy, xx)-s, 1); h{yy, xx}'; zeros(abs(length(tnew) - lidx(yy, xx)), 1)];
            end
        else
            Htemp{yy, xx}  = [zeros(fidx(yy, xx)-1, 1); h{yy, xx}'; zeros(abs(length(tnew) - lidx(yy, xx)), 1)];
        end
		Htotal = Htotal + Htemp{yy, xx};
    end
end

figure(2)
for index = 1:Nscanned
plot(t{index}, h{index})
pause(10/Nscanned)
hold on
end

figure(3)
plot(tnew, Htotal)
hold off