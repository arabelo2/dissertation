tic;
% profile on -historysize 100000000
% clear all
% close all

% % Adds the specified folders to the top of the search path for the current MATLAB® session.
addpath('E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\code\src\Array\', ...
        'E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\code\src\Rectangular\', ...
        'E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\code\src\Circular\UF-Program\', ...
        'E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\code\src\Circular\', ...
        'E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Temp\PPGEM\Dissertação\Programa\Matuda\', ...
        'D:\backup\Usuários\rabeloal\Program\Experimento\202405\Dados_Ensaio_1MHz_20_05_2024\Dados_Ensaio_1MHz_20_05_2024\1MHz_pt_1500\Diretividade')

CW = 0; % CW = 1 : Excitação harmônica (0 -- > Desabilita)
phase = 0; % Phase shift
D = 0.01905; % Diameter [m]
R = D/2; % Radius [m]
c1 = 1500; % [m/s]
rho = 1000; % Density of liquid water [m3/kg]
f0 = 1e6; % Operating frequency of the circular transducer [Hz]
fs = 1.920030720491530e+07; % Sample frequency [Hz] | r 1249
% fs = 9.600038400153603e+06; % Sample frequency [Hz] | az 1500
lambda = c1 / f0; % [m]
STEP = lambda / 10;
NF = R^2/lambda;% Near Field Length or Transition from Near Field to Far Field
k = 2*pi/lambda; % wave number
Uo = 1; % [V]
K = 1; % Constant of the output voltage
nc = 1; % Number of cycles
Dhydrophone = 1e-3; % | r 1249;
% Dhydrophone = 1.5e-3; % | az 1500
discretization = 1;
yc = 0;

Dh = Dhydrophone;
Dt = D;

xmin = 0;  % z-axis
% xmax = 0.210; % z-axis % | r 1249
% xmax = 0.180; % z-axis % | az 1500

ymin = -0.030; % x-axis
ymax = +0.030; % x-axis

xpoints = 71; % | r 1249
% xpoints = 61;  % | az 1500

ypoints = 13;

% x = linspace(xmin, xmax, xpoints); x = x + 4.8375e-04;
x = [0 10 20 30 50 60 80 100 120]*1e-3; x = x + 4.8375e-04;
y = linspace(ymin, ymax, ypoints);

% Velocity potential impulse response of rectangular pistonlike transducers
h_c = cell(length(y), length(x));
t_c = cell(length(y), length(x));

% Piston velocity excitation pulses        
% The input velocity of disc surface
% t_vn = 0 : 1/fs : nc*1/f0;
% vn = InputVelocityOfDiscSurface (Uo, f0, t_vn, phase);
% vn = X_estimated';

% N = 13; [~, V, ~] = q_sequence(N);

% % 1249
% cd 'E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\Experimento\202405\Dados_Ensaio_1MHz_20_05_2024\Dados_Ensaio_1MHz_20_05_2024\1MHz_pt_1249\Perfil_Acustico_x0_z2_z210_Passo3mm'
% d=dir(fullfile((pwd),'*.csv'));
% N = length(d);
% [~, V, ~] = r_sequence(N);

% cd 'E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\Experimento\202405\Dados_Ensaio_1MHz_20_05_2024\Dados_Ensaio_1MHz_20_05_2024\1MHz_pt_1249\Directividade'
N = 13;
[~, V, ~] = d_sequence(N);
% Vvmax = max(max(Vv));

% % 1500
% cd 'E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\Experimento\202405\Dados_Ensaio_1MHz_20_05_2024\Dados_Ensaio_1MHz_20_05_2024\1MHz_pt_1500\Perfil_Acustico_x0_z2_120mm_Passo3mm'
% d=dir(fullfile((pwd),'*.csv'));
% N = length(d);
% [t, V, Vv] = az_sequence(N);

% cd 'E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\Experimento\202405\Dados_Ensaio_1MHz_20_05_2024\Dados_Ensaio_1MHz_20_05_2024\1MHz_pt_1500\Diretividade'
% N = 13;
% [~, V, ~] = q_sequence(N);
% Vvmax = max(max(Vv));

for xx = 1:length(x)
    for yy = 1:length(y)
        xc = y(yy);        
        zc = x(xx);
        [t_temp, h_temp] = summationc(xc, yc, zc, f0, fs, c1, Dh, Dt, discretization);
        h_c{yy, xx} = h_temp;
        t_c{yy, xx} = t_temp;
        
        %%%% texcitation{xx, yy, zz} = texcitation_temp;
     
         y_output = V{yy, xx};
         Y_output = fft(y_output);
         H = fft([h_temp; zeros(length(y_output) - length(h_temp), 1)]);
         X_estimated = ifft(Y_output ./ H);
        
         % v_c{yy, xx} = vn;
         % v_c{yy, xx} = X_estimated;        

         % Digital Filtering: Low-Pass Filter
         % Design a low-pass filter to remove high-frequency noise from the signal

         % Create a low-pass filter (adjust parameters as needed)

         cutoff_freq = 50; % Cutoff frequency (Hz)
         filter_order = 6; % Filter order
         lp_filter = designfilt('lowpassfir', 'FilterOrder', filter_order, 'CutoffFrequency', cutoff_freq, 'SampleRate', fs);

         % Apply the filter to the noisy signal
         filtered_signal = filter(lp_filter, X_estimated);
         v_c{yy, xx} = filtered_signal;

%          % Plot the original and filtered signals
%          figure(70);
%          t_orig = 1/fs*(0:length(X_estimated) - 1);
%          plot(t_orig, X_estimated/max(X_estimated), 'b', 'LineWidth', 1.5);
%          ylim([-1 1])
%          hold on;
%          plot(t_orig, filtered_signal/max(filtered_signal), 'r', 'LineWidth', 2);
%          xlabel('Time');
%          ylabel('Signal');
%          title('Estimated Signal vs. Filtered Signal');
%          legend('Estimated Signal', 'Filtered Signal');
%          grid on;
%          grid minor;
%          hold off;         

%          figure(80)
%          t_e = 1/fs*(0:length(y_output) - 1);
%          plot(t_e, y_output/max(y_output))
%          ylim([-1 1])
%          grid on;
%          grid minor;
%          hold off;
%          pause(.1)

        % Pressure distribution
        Pi_temp = rho*diff(h_temp)/(t_temp(2) - t_temp(1)); 
        % Pi_temp = rho*diff(h_temp); 
        Pi_c{yy, xx} = Pi_temp;
            
        % Pressure waveform
        % Transient pressure
        % p_temp = rho*conv(h_temp, diff(vn)/(t_temp(2) - t_temp(1)));
        
        % p_temp = rho*conv(h_temp, diff(vn));
        p_temp = rho*conv(h_temp, diff(v_c{yy, xx}));
        
        % p_temp = rho*conv(h_temp, vn);        
        % p_temp = conv(vn, Pi_temp);

        P_c{yy, xx} = p_temp;
        t_conv_temp = t_temp(1) + (t_temp(2) - t_temp(1))*(0:1:length(p_temp)-1);
        t_conv_c{yy, xx} = t_conv_temp;
               
        if CW % Artifício para calcular o campo em excitação contínua usando burst longo
            
            %[TF, S1, S2] = ischange(p_temp, "linear");
            %         mid = floor(length(S1)/2);
            %         findex = findchangepts(diff(S1(1:mid)),'Statistic','rms');
            %         lindex = mid + findchangepts(diff(S1(mid+1:end)),'Statistic','rms');
            %         p_temp(1:findex) = 0;
            %         p_temp(lindex:end) = 0;
            
            N = round(fs/f0); % Sample
            mid = floor(length(p_temp)/2);

            if mid > N
                index = mid - N : 1 : mid + N; % garantindo 2 ciclos para buscar max(abs(p))
                p_temp = p_temp(index); % vetor das pressões somente da parte central (2 ciclos) da onda completa
            end
        end

        % Peak amplitude
        Pp_c(yy, xx) = max(abs(p_temp));
            
        % Peak-to-peak amplitude
        % Pp_c(yy, xx) = abs(max(p_temp) - min(p_temp)); 
            
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
% grid on
% grid minor
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)

% figure(30)
% pcolor(z, x, Pp_c/max(max(Pp_c)))
% xlabel('z(m)', 'Color', 'k', 'interpreter', 'latex')
% ylabel('x(m)', 'Color', 'k', 'interpreter', 'latex')
% title(['Pressão Normalizada (Método Analítico) com', num2str(nc),' ciclo(s)'], 'Color', 'k', 'interpreter', 'latex')
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

% figure(31)
% mesh(z, x, Pp_c/max(max(Pp_c)))
% xlabel('z(m)', 'Color', 'k', 'interpreter', 'latex')
% ylabel('x(m)', 'Color', 'k', 'interpreter', 'latex')
% zlabel('Pressão Normalizada', 'Color', 'k', 'interpreter', 'latex')
% title(['Método Analítico com', num2str(nc),' ciclo(s)'], 'Color', 'k', 'interpreter', 'latex')
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

% figure(10)
% hold on
% plot(z*lambda/(R^2), Pp_c_prg(floor(length(x)/2) + 1, :)/max(max(Pp_c_prg)))
% ylabel('Pressão normalizada', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% xlabel('$$z \frac{\lambda}{a^2}$$', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% title('Pressão ao longo do eixo acústico (eixo z)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% grid on
% grid minor

% figure(2)
% hold on
% plot(z*lambda/(R^2), Pp_c_prg(floor(length(x)/2) + 1, :)/max(max(Pp_c_prg)), 'r')
% plot(z*lambda/(R^2), Pp_c_prg(floor(length(x)/2) + 1, :)/max(max(Pp_c_prg)), 'k*')
% ylabel('Pressão normalizada', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% xlabel('$$z \frac{\lambda}{a^2}$$', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% % title('Pressão p/ x=9,5250mm (Raio do transdutor - Dt = 3/4")', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% % title('Pressão ao longo do eixo acústico (eixo z)', 'FontSize', 20, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% title({
% ['Pressão ao longo do eixo acústico (eixo z) para fo = ' num2str(f0/1e6) ' MHz, fs = ' num2str(fs/1e6) ' MHz, x = 0 mm e diâmetro do hidrofone = ' num2str(Dhydrophone*1000) ' mm' ]
% ['com discretização através de uma malha ' num2str(discretization) 'x' num2str(discretization) ' elementos para o levantamento analítico']
% });
% set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
% grid on
% grid minor
% hold off

for index = 1:size(z,2)
figure(index)
hold on
plot(x, Pp_c_prg(:, index)/max(max(Pp_c_prg)), "r")
plot(x, Pp_c_prg(:, index)/max(max(Pp_c_prg)), 'k*')
ylabel('Pressão normalizada', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
xlabel('x (mm)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% title('Excitação Harmônica - perfil de pressão p/ z=90mm (pos. 1)')
title({
['Perfil de pressão acústica (eixo x) para z = ' num2str(z(index)*1000) ' mm, fs = ' num2str(fs/1e6) ' MHz e diâmetro do hidrofone = ' num2str(Dhydrophone*1000) ' mm' ]
['com discretização através de uma malha ' num2str(discretization) 'x' num2str(discretization) ' elementos para o levantamento analítico']
});
legend('Experimental: r 1249', '','Analítico: r 1249','')
set(gca, 'FontName', 'Times New Roman', 'FontSize', 20)
grid minor
grid on
end

% figure()
% plot(x*1000, Pp_c_prg(:, 387)/max(max(Pp_c_prg)))
% ylabel('Pressão normalizada', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% xlabel('x (mm)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% title('Excitação Harmônica - perfil de pressão p/ z=~60 mm (pos. 387)')
% grid on
% grid minor


% figure()
% plot(x*1000, Pp_c_prg(:, 521)/max(max(Pp_c_prg)))
% ylabel('Pressão normalizada', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% xlabel('x (mm)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% title('Excitação Harmônica - perfil de pressão p/ z=~80 mm (pos. 521)')
% grid on
% grid minor

% figure()
% plot(x*1000, Pp_c_prg(:, 654)/max(max(Pp_c_prg)))
% ylabel('Pressão normalizada', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% xlabel('x (mm)', 'FontSize', 14, 'FontWeight', 'bold', 'Color', 'k', 'interpreter', 'latex')
% title('Excitação Harmônica - perfil de pressão p/ z=~100 mm (pos. 654)')
% grid on
% grid minor

% figure()
% plot(z*lambda/(R^2), Pp_c_prg(floor(length(x)/2) + 1, :)/max(max(Pp_c_prg)))

% profile viewer 
toc;