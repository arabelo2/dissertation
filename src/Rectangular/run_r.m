tic;
a = 3/1000;
b = 9.5/1000;
c = 1480;
c1= c;
rho = 1000;
f0 = 2.25e6; % Transducer center frequency [Hz]
fs=100e6; % Sampling frequency [Hz]
lambda = c / f0;
STEP = lambda/3;
NF = a^2/lambda;% Near Field Length or Transition from Near Field to Far Field
x = [-.015:STEP:.015];
y = 0;
z = [0:STEP:3*NF];

% Velocity Potential Impulse Response Of Rectangular Pistonlike Transducers
h_r = cell(length(y), length(x), length(z));
t_r = cell(length(y), length(x), length(z));
%%%%%texcitation = cell(length(x), length(y), length(z));
v = cell(length(y), length(x), length(z));
Pi = cell(length(y), length(x), length(z));
P = cell(length(y), length(x), length(z));
t_conv = cell(length(y), length(x), length(z));
Pp = zeros(length(y), length(x));
Ppp = zeros(length(y), length(x));
Prms = zeros(length(y), length(x));

for zz = 1:length(z)
    for yy = 1:length(y)      
        for xx = 1:length(x)
        [t_temp, h_temp] = vpirOfRectangularPistonlikeTransducers(a, b, c, x(xx), y(yy), z(zz), fs);
        h_r{yy, xx, zz} = h_temp;
        t_r{yy, xx, zz} = t_temp;
        
        % Piston velocity excitation pulses
        % Excitation - Sitau + array of 5MHz - Format: HDF5
        v_temp = resample(hdf5read('ref_pulse-40MHz.h5', 'ascan'), fs*1e-6, 40);
        
        % FILTER
        v_temp = v_temp.*hanning(max(size(v_temp)));
        C = 1/max(v_temp);
        v_temp = C*v_temp;   
        
        %%%% texcitation{xx, yy, zz} = texcitation_temp;
        v{yy, xx, zz} = v_temp;
            
        % Pressure distribution
        Pi_temp = rho*diff(h_temp)/(t_temp(2) - t_temp(1)); 
        Pi{yy, xx, zz} = Pi_temp;
            
        % Transient pressure
        p_temp = rho*conv(h_temp, diff(v_temp)/(t_temp(2) - t_temp(1)));
        % p_temp = conv(v_temp, Pi_temp);
        P{yy, xx, zz} = p_temp;
        t_conv_temp = t_temp(1) + (t_temp(2) - t_temp(1))*(0:1:length(p_temp)-1);
        t_conv{yy, xx, zz} = t_conv_temp;
            
        % Peak amplitude
        Pp(xx, zz) = max(p_temp);
            
        % Peak-to-peak amplitude
        Ppp(xx, zz) = max(p_temp) - min(p_temp); 
            
        % Root mean square amplitude
        Prms(xx, zz) = (max(p_temp) - min(p_temp))/(2*sqrt(2));
        end
    end
end

figure()
pcolor(z, x, Ppp)
shading interp
colormap(jet)
colorbar
axis normal
ylabel('x(mm)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
xlabel('z(mm)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
title('Pressure field',  'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
grid on
grid minor
set(gca,'Ydir','reverse');

figure()
mesh(z, x, Ppp)
shading interp
colormap(jet)
colorbar
axis normal
ylabel('x(mm)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
xlabel('z(mm)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
title('Pressure field',  'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
grid on
grid minor
set(gca,'FontSize',20);
toc;
