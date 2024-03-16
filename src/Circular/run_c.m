tic;
profile on -historysize 100000000
%clear all
%close all

% % Adds the specified folders to the top of the search path for the current MATLAB® session.
addpath('E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\code\src\Array\', ...
        'E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\code\src\Rectangular\', ...
        'E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\code\src\Circular\UF-Program\', ...
        'E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\code\src\Circular\', ...
        'E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Temp\PPGEM\Dissertação\Programa\Matuda\')

D = 0.01905; % Diameter [m]
R = D/2; % Radius [m]
c1 = 1500; % [m/s]
rho = 1000; % Density of liquid water [m3/kg]
f0 = 1e6; % Operating frequency of the circular transducer [Hz]
fs = 32e6; % Sample frequency [Hz]
lambda = c1/f0; % [m]
STEP = lambda/10;
NF = R^2/lambda;% Near Field Length or Transition from Near Field to Far Field
k = 2*pi/lambda; % wave number
Uo = 1; % [V]
K = 1; % Constant of the output voltage
nc = 27; % Number of cycles
%Dhydrophone = 0.6e-3;
%discretization = 7;

xmin = 0.080;  % z-axis
xmax = 0.100; % z-axis

ymin = -0.030; % x-axis
ymax = 0.030; % x-axis

xpoints = 2;
ypoints = 407;

x = linspace(xmin, xmax, xpoints);
y = linspace(ymin, ymax, ypoints);

% Velocity potential impulse response of rectangular pistonlike transducers
h_c = cell(length(y), length(x));
t_c = cell(length(y), length(x));

% Piston velocity excitation pulses        
% The input velocity of disc surface
t_vn = 0 : 1/fs : nc*1/f0;
vn = InputVelocityOfDiscSurface (Uo, f0, t_vn);

for xx = 1:length(x)
    for yy = 1:length(y)
        [t_temp, h_temp] = vpirOfCircularPistonlikeTransducers(x(xx), y(yy), D, c1, fs);
        h_c{yy, xx} = h_temp;
        t_c{yy, xx} = t_temp;
        
        %%%% texcitation{xx, yy, zz} = texcitation_temp;
        v_c{yy, xx} = vn;
            
        % Pressure distribution
        Pi_temp = rho*diff(h_temp)/(t_temp(2) - t_temp(1)); 
        % Pi_temp = rho*diff(h_temp); 
        Pi_c{yy, xx} = Pi_temp;
            
        % Pressure waveform
        % Transient pressure
        % p_temp = rho*conv(h_temp, diff(vn)/(t_temp(2) - t_temp(1)));
        p_temp = rho*conv(h_temp, diff(vn));
        % p_temp = rho*conv(h_temp, vn);        
        % p_temp = conv(vn, Pi_temp);

        P_c{yy, xx} = p_temp;
        t_conv_temp = t_temp(1) + (t_temp(2) - t_temp(1))*(0:1:length(p_temp)-1);
        t_conv_c{yy, xx} = t_conv_temp;
               
%         [TF, S1, S2] = ischange(p_temp, "linear");
%         mid = floor(length(S1)/2);
%         findex = findchangepts(diff(S1(1:mid)),'Statistic','rms');
%         lindex = mid + findchangepts(diff(S1(mid+1:end)),'Statistic','rms');
%         p_temp(1:findex) = 0;
%         p_temp(lindex:end) = 0;

        N = round(fs/f0); % Sample
        mid = floor(length(p_temp)/2);
        index = mid - N : 1 : mid + N; % garantindo 2 ciclos para buscar max(abs(p))
        p_temp = p_temp(index); % vetor das pressões somente da parte central (2 ciclos) da onda completa
                    
        % Peak amplitude
        Pp_c(yy, xx) = max(abs(p_temp));
            
        % Peak-to-peak amplitude
        % Ppp_c(yy, xx) = abs(max(p_temp) - min(p_temp)); 
            
        % Root mean square amplitude
        % Prms_c(yy, xx) = abs(max(p_temp) - min(p_temp))/(2*sqrt(2));
    end
end

% Axial pressure amplitude for a baffled circular plane piston.
OnAxialPressure = AxialPressure(k, x, R);

z = x;
x = y;

% figure()
% plot(z*lambda/(R^2), OnAxialPressure)
% title(['Axial pressure amplitude for a baffled circular plane piston with ', num2str(nc),' cycle(s)'], 'Color', 'k', 'interpreter', 'latex')
% xlabel('z (mm)')
% ylabel('Normalized pressure','interpreter','latex')
% %ylabel('$\frac{P(r, t)}{2 \rho c U_o}$','interpreter','latex')
% grid on
% grid minor
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)

% figure()
% pcolor(z, x, Pp_c/max(max(Pp_c)))
% xlabel('z(m)', 'Color', 'k', 'interpreter', 'latex')
% ylabel('x(m)', 'Color', 'k', 'interpreter', 'latex')
% title(['Pressure field of a circular plane rigid baffled piston with ', num2str(nc),' cycle(s)'], 'Color', 'k', 'interpreter', 'latex')
% az = 0; % az = -90/90 -- > Horizontal ; % az = 0/180 -- > Vertical;
% el = 90;
% view(az, el);
% shading interp
% colormap(jet)
% colorbar
% axis padded
% grid on
% grid minor
% set(gca,'Ydir','reverse');
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
% daspect('auto') % daspect([1 1 1])


% figure()
% mesh(z, x, Pp_c/max(max(Pp_c)))
% xlabel('z(m)', 'Color', 'k', 'interpreter', 'latex')
% ylabel('x(m)', 'Color', 'k', 'interpreter', 'latex')
% % zlabel('Max(|P(r,t)|)')
% zlabel('Normalized pressure', 'Color', 'k', 'interpreter', 'latex')
% title(['Pressure field of a circular plane rigid baffled piston with ', num2str(nc),' cycle(s)'], 'Color', 'k', 'interpreter', 'latex')
% shading interp
% colormap(jet)
% colorbar
% axis padded
% grid on
% grid minor
% set(gca,'Ydir','reverse');
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
% daspect('auto') % daspect([1 1 1])

Pp_c_prg = Pp_c;

% figure()
%plot(z*lambda/(R^2), Pp_c_prg(floor(length(x)/2) + 1, :)/max(max(Pp_c_prg)))
% ylabel('Pressão normalizada', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% xlabel('$$z \frac{\lambda}{a^2}$$', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% title('Pressão ao longo do eixo acústico (eixo z)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% grid on
% grid minor

%figure(4)
%hold on
%plot(x*1000, Pp_c_prg(:, 1)/max(max(Pp_c_prg)))

% figure(10)
% plot(x*1000, Pp_c_prg(:, 121)/max(max(Pp_c_prg)))
% figure(11)
% plot(x*1000, Pp_c_prg(:, 387 )/max(max(Pp_c_prg)))
% figure(12)
% plot(z*lambda/(R^2), Pp_c_prg(floor(length(x)/2) + 1, :)/max(max(Pp_c_prg)))
profile viewer 
toc;
