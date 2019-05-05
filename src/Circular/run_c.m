tic;
D = 19/1000;
c = 1480;
c1= c;
rho = 1000;
f0 = 2.25e6; % Transducer center frequency [Hz]
fs = 100e6; % Sampling frequency [Hz]
lambda = c / f0;
STEP = lambda/3;
NF = D^2/4/lambda;% Near Field Length or Transition from Near Field to Far Field
x = [0:STEP:3*NF];
y = [-.015:STEP:.015];

% Velocity potential impulse response of rectangular pistonlike transducers
h_c = cell(length(y), length(x));
t_c = cell(length(y), length(x));

for xx = 1:length(x)
    for yy = 1:length(y)
        [t_temp, h_temp] = vpirOfCircularPistonlikeTransducers(x(xx), y(yy), D, c, fs);
        h_c{yy, xx} = h_temp;
        t_c{yy, xx} = t_temp;
        
        % Piston velocity excitation pulses
        % Excitation - Sitau + array of 5MHz - Format: HDF5
        v_temp = resample(hdf5read('ref_pulse-40MHz.h5', 'ascan'), fs*1e-6, 40);
        
        % FILTER
        v_temp = v_temp.*hanning(max(size(v_temp)));
        C = 1/max(v_temp);
        v_temp = C*v_temp;   
        
        %%%% texcitation{xx, yy, zz} = texcitation_temp;
        v{yy, xx} = v_temp;
            
        % Pressure distribution
        Pi_temp = rho*diff(h_temp)/(t_temp(2) - t_temp(1)); 
        Pi{yy, xx} = Pi_temp;
            
        % Transient pressure
        p_temp = rho*conv(h_temp, diff(v_temp)/(t_temp(2) - t_temp(1)));
        % p_temp = conv(v_temp, Pi_temp);
        P{yy, xx} = p_temp;
        t_conv_temp = t_temp(1) + (t_temp(2) - t_temp(1))*(0:1:length(p_temp)-1);
        t_conv{yy, xx} = t_conv_temp;
            
        % Peak amplitude
        Pp(yy, xx) = max(p_temp);
            
        % Peak-to-peak amplitude
        Ppp(yy, xx) = max(p_temp) - min(p_temp); 
            
        % Root mean square amplitude
        Prms(yy, xx) = (max(p_temp) - min(p_temp))/(2*sqrt(2));
    end
end

figure()
pcolor(x, y, Ppp)
shading interp
colormap(jet)
colorbar
axis normal
ylabel('x(mm)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
xlabel('y(mm)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
title('Pressure field',  'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
grid on
grid minor
set(gca,'Ydir','reverse');

figure()
mesh(x, y, Ppp)
shading interp
colormap(jet)
colorbar
axis normal
ylabel('y(mm)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
xlabel('x(mm)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
title('Pressure field',  'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
grid on
grid minor
set(gca,'FontSize',20);
toc;
