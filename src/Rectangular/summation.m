a = 6e-3;
b = 9.5e-3;
c1 = 1500;
fs = 64e6;
xc = 0e-3;
yc= 0e-3;
zc = 1e-3;
D = 1e-3;
discretization = 21;

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
grid minor
grid on
hold on
pause(60/Nscanned)
grid off
end

figure(2)
plot(tnew, Htotal)
grid minor
grid on
hold off
grid off
