tic;
clear all
close all
D = 0.01905; % Diameter [m]
R = D/2; % Radius [m]
c1 = 1500; % [m/s]
rho = 1000; % Density of liquid water [m3/kg]
f0 = 1e6; % Operating frequency of the circular transducer [Hz]
fs = 32e6; % Sample frequency [Hz]
lambda = c1 / f0; % [m]
STEP = lambda / 5;
NF = D^2/4/lambda;% Near Field Length or Transition from Near Field to Far Field
k = 2*pi*lambda^-1; % wave number
Uo = 1; % [V]
K = 1; % Constant of the output voltage
nc = 17; % Number of cycles
Dhydrophone = 0.0002;
discretization = 9;

xmin = 0.002;  % z-axis
xmax = +0.202; % z-axis

ymin = -0.020; % x-axis
ymax = +0.020; % x-axis

xpoints = 401; 
ypoints = 81;

x = linspace(xmin, xmax, xpoints);
y = linspace(ymin, ymax, ypoints);

% Velocity potential impulse response of rectangular pistonlike transducers
h_c = cell(length(y), length(x));
t_c = cell(length(y), length(x));

for xx = 1:length(x)
    for yy = 1:length(y)
        [t_temp, h_temp] = vpirOfCircularPistonlikeTransducersWScanner(x(xx), y(yy), D, c1, fs, Dhydrophone, discretization);
        h_c{yy, xx} = h_temp;
        t_c{yy, xx} = t_temp;
        
        % Piston velocity excitation pulses        
        % The input velocity of disc surface
        t_vn = 0 : 1/fs : nc*1/f0;
        vn = InputVelocityOfDiscSurface (Uo, f0, t_vn);

        % FILTER
        % vn = vn.*hanning(max(size(vn)));
        % C = 1/max(vn);
        % vn = C*vn;   
        
        %%%% texcitation{xx, yy, zz} = texcitation_temp;
        v_c{yy, xx} = vn;
            
        % Pressure distribution
        Pi_temp = rho*diff(h_temp)/(t_temp(2) - t_temp(1)); 
        Pi_c{yy, xx} = Pi_temp;
            
        % Transient pressure
        p_temp = rho*conv(h_temp, diff(vn)/(t_temp(2) - t_temp(1)));
        
        % p_temp = conv(vn, Pi_temp);
        P_c{yy, xx} = p_temp;
        t_conv_temp = t_temp(1) + (t_temp(2) - t_temp(1))*(0:1:length(p_temp)-1);
        t_conv_c{yy, xx} = t_conv_temp;
            
        % Peak amplitude
        Pp_c(yy, xx) = max(abs(p_temp));
            
        % Peak-to-peak amplitude
        Ppp_c(yy, xx) = abs(max(p_temp) - min(p_temp)); 
            
        % Root mean square amplitude
        Prms_c(yy, xx) = abs(max(p_temp) - min(p_temp))/(2*sqrt(2));
    end
end

% Axial pressure amplitude for a baffled circular plane piston.
OnAxialPressure = AxialPressure(k, x, R);

z = x;
x = y;

figure()
plot(z, OnAxialPressure)
title(['Axial pressure amplitude for a baffled circular plane piston with ', num2str(nc),' cycle(s)'], 'Color', 'k', 'interpreter', 'latex')
xlabel('z (mm)')
ylabel('Normalized pressure','interpreter','latex')
%ylabel('$\frac{P(r, t)}{2 \rho c U_o}$','interpreter','latex')
grid on
grid minor
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)

figure()
pcolor(z, x, Pp_c/max(max(Pp_c)))
xlabel('z(mm)', 'Color', 'k', 'interpreter', 'latex')
ylabel('x(mm)', 'Color', 'k', 'interpreter', 'latex')
title(['Pressure field of a circular plane rigid baffled piston with ', num2str(nc),' cycle(s)'], 'Color', 'k', 'interpreter', 'latex')
az = 0; % az = -90/90 -- > Horizontal ; % az = 0/180 -- > Vertical;
el = 90;
view(az, el);
shading interp
colormap(jet)
colorbar
axis padded
grid on
grid minor
set(gca,'Ydir','reverse');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
daspect('auto') % daspect([1 1 1])

figure()
mesh(z, x, Pp_c/max(max(Pp_c)))
xlabel('z(mm)', 'Color', 'k', 'interpreter', 'latex')
ylabel('x(mm)', 'Color', 'k', 'interpreter', 'latex')
% zlabel('Max(|P(r,t)|)')
zlabel('Normalized pressure', 'Color', 'k', 'interpreter', 'latex')
title(['Pressure field of a circular plane rigid baffled piston with ', num2str(nc),' cycle(s)'], 'Color', 'k', 'interpreter', 'latex')
shading interp
colormap(jet)
colorbar
axis padded
grid on
grid minor
set(gca,'Ydir','reverse');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
daspect('auto') % daspect([1 1 1])

toc;
