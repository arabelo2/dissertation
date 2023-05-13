%%%% 
%%%% Loading data
%%%% 
%%%% Required input
%
% a -- > Half length in [m].
% b -- > Half width/height in [m].
% c -- > The velocity of sound in the propagating medium [m/s].
% x -- > Abscissa (X-axis) [m].
% y -- > Ordinate (Y-axis) [m].
% z -- > Applicate (Z-axis) [m] -- > The Z-axis is perpendicular to the plane XY.
% f0 -- > The operating frequency of US-transducer [Hz]
% fs -- > The sample frequency [Hz]

% Adds the specified folders to the top of the search path for the current MATLAB® session.
addpath('E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\code\src\Array')
addpath('E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\code\src\Rectangular')
addpath('E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Temp\PPGEM\Dissertação\Programa\Matuda')

% Start stopwatch timer
tic

% Remove items from workspace
% clear all;
% Delete all figures
close all;

format longG;
%format compact;

% Data
rho = 1000; % [kg/m^3]
c1 = 1500; % [m/s]
f0 = 3e6; % [Hz]
f = f0;
fs=16*f0; % Sampling frequency [Hz]
lambda = c1 / f0; % [m]
kerf = 6e-4; % [m]
STEP = lambda/3;  % [m]

% X-axis
M = 1; % Number of elements (Columns)
a = 15e-3/2; % Half of width of element [m]

% NF = a^2/lambda;% Near Field Length or Transition from Near Field to Far Field

% Y-axis
N = 1; % Number of elements (Rows)
b = 1.6*a; % m

% Delay law 
delayLawEnabled = 0; % 0 --> OFF and 1 --> ON

xmin = -0.012; % xmin = -(2*a+kerf)*(M/2+1);
xmax = +0.012; % xmax = (2*a+kerf)*(M/2+1);
ymin = 0;
ymax = +0;
zmin = 0.002;
zmax = +0.402; % m -- > The Z-axis is perpendicular to the plane XY.

xpoints = 151;
ypoints = 1;
zpoints = 2401;

x = linspace(xmin, xmax, xpoints);
y = linspace(ymin, ymax, ypoints);
z = linspace(zmin, zmax, zpoints);

z_idx = (z >= .020 & z <= .400 );
x_idx = (x >= -.012 & x <= .012 );

% Focal distance
F = 10000; % [m]

% Rotate around the z-axis (Roll)
PHII = 0; % Degree angle

% az - Rotate around the y-axis (Pitch) / Azimuth is the counterclockwise angle in the z-x plane measured in radians from the positive z-axis.
THETA = 0; % Degree angle

% el - Rotate around the x-axis (Yaw) / Elevation is the elevation angle in radians from the z-x plane 
PSI = 0; % Degree angle

% Transforms corresponding elements of the focal point spherical coordinate arrays azimuth (theta), elevation (psi), and F to Cartesian.
[zf, xf, yf] = focalPoint (THETA, PSI, F);
 
% Velocity Potential Impulse Response Of Rectangular Pistonlike Transducers
h = cell(length(y), length(x), length(z));
t = cell(length(y), length(x), length(z));
%%%%%texcitation = cell(length(x), length(y), length(z));
v = cell(length(y), length(x), length(z));
Pi = cell(length(y), length(x), length(z));
P = cell(length(y), length(x), length(z));
t_conv = cell(length(y), length(x), length(z));
Pp = zeros(length(y), length(x));
Ppp = zeros(length(y), length(x));
Prms = zeros(length(y), length(x));
for zz = 1:length(z(1, :))
    for yy = 1:length(y(1, :))
        for xx = 1:length(x(1, :))
            %[t_temp, h_temp] = vpirOfRectangularPistonlikeTransducers(a, b, c1, x(xx), y(yy), z(zz), f, sample);
            [h_temp, t_temp, td, ex, ey, ez, dDtmn, exm, eyn, B2x, B2y] = vpirOfRectangularArrayPistonlikeTransducers(a, b, c1, x(xx), y(yy), z(zz), fs, N, M, kerf, kerf, delayLawEnabled, zf, xf, yf, F);
            h{yy, xx, zz} = h_temp;
            t{yy, xx, zz} = t_temp;

            % Piston velocity excitation pulses
            % Wideband, type I pulse.
%             K = 3.833;
%             C = 1; 
%             texcitation_temp = 0;
%             count = 0;
%             while c1*texcitation_temp(end) < 0.00300
%                 texcitation_temp(count + 1) = count*(t_temp(2) - t_temp(1));
%                 count = count + 1; 
%             end            
            % Narrow-band, type II pulse.
            K = 1.437;
            C = 1;
            texcitation_temp = 0;
            count = 0;
            while c1*texcitation_temp(end) < 0.00600
                texcitation_temp(count + 1) = count*(t_temp(2) - t_temp(1)); 
                count = count + 1;
            end

            v_temp = C*texcitation_temp.^3.*exp(-K*f*texcitation_temp).*cos(2*pi*f*texcitation_temp);             
            % CW Excitation (But it's not working.)
            % v_temp = cos(2*pi*f*texcitation_temp);                                  
            texcitation{xx, yy, zz} = texcitation_temp;
            v{xx, yy, zz} = v_temp;
            
            % Pressure distribution
            Pi_temp = rho*diff(h_temp)/(t_temp(2) - t_temp(1)); 
            Pi{yy, xx, zz} = Pi_temp;
            
            % Transient pressure
            p_temp = rho*conv(h_temp, diff(v_temp)/(t_temp(2) - t_temp(1)), 'full');
            % p_temp = conv(v_temp, Pi_temp);
            P{yy, xx, zz} = p_temp;
            t_conv_temp = t_temp(1) + (t_temp(2) - t_temp(1))*(0:1:length(p_temp)-1);
            t_conv{yy, xx, zz} = t_conv_temp;
            
            % Peak amplitude
            Pp(xx, zz) = max(abs(p_temp));
            
            % Peak-to-peak amplitude
            Ppp(xx, zz) = max(p_temp) - min(p_temp); 
            
            % Root mean square amplitude
            Prms(xx, zz) = (max(abs(p_temp)) - min(abs(p_temp)))/(2*sqrt(2));            
        end
    end
end

% - %%%%%%%%%%%%%%%%%%%%%
%   PLOT
% - %%%%%%%%%%%%%%%%%%%%%

% figure(2)
% subplot(1,2,1)
% plot(exm, eyn, 'r+', ... % The centroid of the mth element
%     exm + a, eyn + b, 'b.', ...
%     exm + 0, eyn + b, 'b.', ...
%     exm + a, eyn + 0, 'b.', ...
%     exm - a, eyn - b, 'b.', ...
%     exm + 0, eyn - b, 'b.', ...
%     exm - a, eyn + 0, 'b.', ...
%     exm - a, eyn + b, 'b.', ...
%     exm + a, eyn - b, 'b.', ...
%     x, y, 'k*')
% axis equal
% xlabel('X-axis (m)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% ylabel('Y-axis (m)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% title('Geometrical region (Array)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% grid on
% grid minor
% set(gca,'FontSize',20);

% figure(2)
% subplot(1,2,2)
% plot(t{1}*c1, h{1}/c1)
% % figure(4)
% % subplot(2,2,1)
% % plot(t{1}*c1, h{1}/c1)
% axis normal
% xlabel('ct', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% %xlabel('ct (m)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% ylabel('$$\frac{h(\overline{x}, t)}{c}$$',  'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% title('Temporal evolution of the impulse response', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% grid on
% grid minor
% set(gca,'FontSize',20);

% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Temporal evolution of the impulse response (mono) - d', 'jpg')
% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Temporal evolution of the impulse response (mono) - d', 'bmp')
% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Temporal evolution of the impulse response (mono) - d', 'fig')

% Plot individual velocity potential impulse response
% for r = 1:length(x(1, :))
%     for s = 1:length(y(1, :))
%         for u = 1:length(z(1, :))
%             plot(c1*t{r, s, u}/a, h{r, s, u}/c1, '.')        
%             axis normal
%             xlabel('ct/a', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
%             ylabel('h($\overline{x}$, t)/c', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
%             title('Temporal evolution of the impulse response of wide rectangular source (Array)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')       
%             grid on
%             grid minor
%             set(gca,'FontSize',20);
%             pause(.25);
%             hold on
%         end
%     end
% end

% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Temporal evolution of the impulse response of narrow rectangular source (mono) - b', 'jpg')
% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Temporal evolution of the impulse response of narrow rectangular source (mono) - b', 'bmp')
% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Temporal evolution of the impulse response of narrow rectangular source (mono) - b', 'fig')


% Piston velocity excitation pulses v(t)
% figure(3)
% plot(c1*texcitation{1}*1000, v{1})
% axis normal
% xlabel('ct(mm)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% ylabel('v(t)(m/s)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% title('Piston velocity excitation pulses - Wideband (Type I)',  'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% % title('Piston velocity excitation pulses - Narrow-band (Type II)',  'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% grid on
% grid minor
% set(gca,'FontSize',20);

% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Piston velocity excitation pulses - Narrow-band (Type II)', 'jpg')
% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Piston velocity excitation pulses - Narrow-band (Type II)', 'bmp')
% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Piston velocity excitation pulses - Narrow-band (Type II)', 'fig')


% %Pressure distribution
% figure(4)
% %subplot(2,2,2)
% %plot(c1*t{1}(1:cellfun(@length, t(1)) - 1)/a, Pi{1}*(t{1}(2) - t{1}(1))/rho/c1)
% plot(c1*t{1}(1:cellfun(@length, t(1)) - 1)*1000, Pi{1}*(t{1}(2) - t{1}(1))/rho/c1)
% % plot(c1*t{1}(1:cellfun(@length, t(1)) - 1)/a, Pi{1}/mean(sqrt(Pi{1}.^2))/4)
% % plot(c1*t{1}(1:cellfun(@length, t(1)) - 1)/a, Pi{1}/mean(sqrt(Pi{1}.^2)))
% axis normal
% %xlabel('ct/a', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% xlabel('ct (mm)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% ylabel('$$\frac{Pi(\overline{x}, t)\partial{t}}{\rho c}$$', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% title('Temporal evolution of the pressure impulse response', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% grid on
% grid minor
% set(gca,'FontSize',20);

% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Temporal evolution of the pressure impulse response at field points in the four geometrical region - d', 'jpg')
% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Temporal evolution of the pressure impulse response at field points in the four geometrical region - d', 'bmp')
% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Temporal evolution of the pressure impulse response at field points in the four geometrical region - d', 'fig')

% Transient pressure
% figure(5)
% %figure(4)
% %subplot(2,2,3)
% plot(t_conv{1}*c1*1000, (t{1}(2) - t{1}(1))*P{1}/c1/rho, '.')
% axis normal
% xlabel('ct(mm)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% ylabel('$$\frac{P(\overline{x}, t)\partial{t}}{\rho c}$$', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% title('Temporal evolution of the transient pressure - Type I', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% grid on
% grid minor
% set(gca,'FontSize',20);

% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Impulse responses and transient pressure waveforms', 'jpg')
% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Impulse responses and transient pressure waveforms', 'bmp')
% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Impulse responses and transient pressure waveforms', 'fig')

% On-axis relative peak amplitude of pressure waveforms 
% On-axis relative peak and RMS amplitudes of the pressure waveforms
% figure(6)
% hold on
% plot(z*1000, cellfun(@max, P(:))*(t{1}(2) - t{1}(1))/c1/rho)
% % plot(z*1000, Ppp(1, :)*((t{1}(2) - t{1}(1))/c1/rho), 'b.')
% % plot(z*1000, Ppp(1, :), 'r.')
% % plot(z*1000, (t{1}(2) - t{1}(1))*Prms(1, :)/c1/rho, 'r')
% axis normal
% xlabel('z(mm)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% ylabel('$$\frac{Pp(\overline{x}, t, z)\partial{t}}{\rho c}$$', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% % ylabel('$$\frac{Ppp(\overline{x}, t, z)\partial{t}}{\rho c}$$', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% title('On-axis relative peak amplitude of the pressure waveforms - Types I and II',  'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% % title('On-axis relative peak and maximum amplitudes of the pressure waveforms - Types I and II',  'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% grid on
% grid minor
% set(gca,'FontSize',20);
 
% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\On-axis relative peak amplitude of the pressure waveforms - Types I and II - Errata', 'jpg')
% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\On-axis relative peak amplitude of the pressure waveforms - Types I and II - Errata', 'bmp')
% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\On-axis relative peak amplitude of the pressure waveforms - Types I and II - Errata', 'fig')

% figure(17)
% pcolor(x*1e3, z*1e3, abs(Pp')*(t{1}(2) - t{1}(1))/c1/rho)
% % pcolor(z*1e3, x*1e3, abs(Pp)*(t{1}(2) - t{1}(1))/c1/rho)
% % pcolor(z*1e3, y*1e3, cellfun(@max, P)*(t{1}(2) - t{1}(1))/c1/rho)
% shading interp
% colormap(jet)
% colorbar
% axis equal
% ylabel('z(mm)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% xlabel('x(mm)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% % xlabel('z(mm)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% % ylabel('x(mm)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% title('Pressure field - Relative peak amplitude - Type I excitation and y = 0 plane.',  'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% grid on
% grid minor
% % set(gca,'FontSize',20);
% set(gca,'Ydir','reverse', 'FontSize',20);

figure(20)
imagesc(x(x_idx)*1e3, z(z_idx)*1e3, abs(Ppp(x_idx, z_idx)')/max(max(abs(Ppp))))
% figure(200)
% imagesc(x*1e3, z*1e3, abs(Pp')/max(max(abs(Pp))))
az = -90; % az = -90/90 -- > Horizontal ; % az = 0/180 -- > Vertical;
el = 90;
view(az, el);
shading interp
colormap(jet)
colorbar
axis padded
ylabel('z(mm)', 'Color', 'k', 'interpreter', 'latex')
xlabel('x(mm)', 'Color', 'k', 'interpreter', 'latex')
title('Pressão Normalizada - Excitação Tipo II e plano y = 0', 'Color', 'k', 'interpreter', 'latex')
grid on
grid minor
set(gca,'Ydir','reverse');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
daspect('auto') % daspect([1 1 1])


% figure(47)
% pz=length(z);
% px=length(x);
% zMin = z(1)*1e3;
% zMax = z(end)*1e3;
% xMin = x(1)*1e3;
% xMax = x(end)*1e3;
% x_axis_theoretical = linspace(xMin, xMax, px);
% z_axis_theoretical = linspace(zMin, zMax, pz);
% z_idxt = (z_axis_theoretical >= 5 & z_axis_theoretical <= 55 );
% x_idxt = (x_axis_theoretical >= -20 & x_axis_theoretical <= 20 );
% imagesc(x_axis_theoretical(x_idxt), z_axis_theoretical(z_idxt), abs(Ppp(x_idxt, z_idxt))')
% %imagesc(x_axis_theoretical, z_axis_theoretical, abs(Ppp)')
% shading interp
% colormap(jet)
% colorbar
% axis normal
% ylabel('z(mm)', 'FontSize', 25, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% xlabel('x(mm)', 'FontSize', 25, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% title('Teórico',  'FontSize', 25, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% grid on
% grid minor
% set(gca,'Ydir','reverse');
% set(gca,'FontSize',25);
% daspect([1 1 1])

% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Pressure field - Relative peak amplitude - Type I excitation and y = 0 plane - Errata', 'jpg')
% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Pressure field - Relative peak amplitude - Type I excitation and y = 0 plane - Errata', 'bmp')
% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Pressure field - Relative peak amplitude - Type I excitation and y = 0 plane - Errata', 'fig')

% Three-dimensional plot of the relative peak amplitude of the pressure waveforms in the near field of a rectangular piston.
figure(8)
mesh(z(z_idx)*1e3, x(x_idx)*1e3, abs(Ppp(x_idx, z_idx))/max(max(abs(Ppp))))
% figure(80)
% mesh(z*1e3, x*1e3, abs(Pp)/max(max(abs(Pp))))
shading interp
colormap(jet)
colorbar('eastoutside')
axis padded
xlabel('z(mm)', 'Color', 'k', 'interpreter', 'latex')
ylabel('x(mm)', 'Color', 'k', 'interpreter', 'latex')
% zlabel('$$\frac{Pp(\overline{x}, y, z, t)\partial{t}}{\rho c}$$', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
zlabel('Pressão Normalizada', 'Color', 'k', 'interpreter', 'latex')
title('Pressão Normalizada - Excitação Tipo II e plano x = 0', 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
grid on
grid minor
set(gca,'Ydir','reverse');
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
daspect ('auto') % daspect([1 1 1])

% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Three-dimensional plot of the relative peak amplitute of the pressure waveforms in the near field of a rectangular piston - Type I and y = 0 plane - Errata', 'jpg')
% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Three-dimensional plot of the relative peak amplitute of the pressure waveforms in the near field of a rectangular piston - Type I and y = 0 plane - Errata', 'bmp')
% saveas(gcf, 'C:\Temp\PPGEM\Dissertação\Qualificação\Figuras\Original pictures\Three-dimensional plot of the relative peak amplitute of the pressure waveforms in the near field of a rectangular piston - Type I and y = 0 plane - Errata', 'fig')

% - %%%%%%%%%%%%%%%%%%%%%
%   UNDER TESTING
% - %%%%%%%%%%%%%%%%%%%%%

% for q = 1:length(z(1, :))
%     plot(t{q}*c1/a, h{q}/c1)
%     pause (.5)
% end

% for q = 1:length(z(1, :))
%     plot(t_conv{q}*c1/a, (t{1}(2) - t{1}(1))*abs(P{q})/c1/rho)
%     pause (.3)
% end

% Ploting the geometrical region
% for r = 1:length(x(1, :))
%     for s = 1:length(y(1, :))
%         for u = 1:length(z(1, :))
%             figure(1)
%             subplot(1,2,2)
%             plot(t{u}*c1/a, h{u}/c1)
%             figure(1)
%             subplot(1,2,1)
%             plot(0,0,'k+', a,b,'b*', 0,b,'b*', a,0,'b*', -a,-b,'b*', 0,-b,'b*', -a,0,'b*', -a,b,'b*', a,-b,'b*', x(r),y(s),'ro')
%             %xlim([-a - abs(x(r)), a + abs(x(r))])
%             %ylim([-b - abs(y(s)), b + abs(y(s))])
%             axis normal
%             xlabel('X-axis (m)')
%             ylabel('Y-axis (m)')
%             title('Geometrical region (Mono)')
%             grid on
%             grid minor 
%             pause (.5)
%         end
%     end
% end

% for zz = 1:length(z(1, :))
%     for yy = 1:length(y(1, :))
%         for xx = 1:length(x(1, :))
%             plot(1000*c1*t_conv{yy, xx, zz}, abs(P{yy, xx, zz})/max(abs(P{yy, xx, zz})))
%             %xlim([0, 500])
%             %ylim([-1, 1])
%             grid on
%             grid minor
%             pause(.3)
%         end
%     end
% end

% %
% % Beam Steering Through a Planar Interface
% %
% figure(20); plot(exm, dDtmn, 'b.')
% xlabel('exm')
% ylabel('Dtdmn')
% figure(21); plot(eyn, dDtmn, 'b.')
% xlabel('eyn')
% ylabel('Dtdmn')
% 
% figure(22); plot3(exm, eyn, dDtmn, 'b.')
% xlabel('exm')
% ylabel('eyn')
% zlabel('Dtdmn')
% grid on
% grid minor

%Plot individual velocity potential impulse response
% for zz = 1:length(z(1, :))
%     for yy = 1:length(y(1, :))
%         for xx = 1:length(x(1, :))
%             figure (5)
%             % plot(c1*t{yy, xx, zz}, h{yy, xx, zz})
%             plot(t{yy, xx, zz}, h{yy, xx, zz})
%             axis square
%             title('Temporal evolution of the impulse response (Array)')             
%             hold on
%             grid on
%             grid minor
%             pause(.3)
%         end
%     end
% end

% Read elapsed time from stopwatch

% Save workspace variables to file
% tic
% filename = sprintf('campoacustico_%s.mat', datestr(now,'yyyymmddTHHMMSS'))
% save(filename,'-v7.3')
% toc

% Read elapsed time from stopwatch
toc

% ALARM

load handel

for alarm=(1:1)
    sound(y,Fs)
    disp(alarm)
    pause (15)
end