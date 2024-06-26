clear all;
close all;

rho = 1000; %Density of liquid water [m3/kg]
c = 1500; %[m/s]
fc=1e6; %Frequencia Central [Hz]
Fface=1e6; %Frequencia na face do transdutor [Hz]
T=fc^-1; %[s]
lambda=c*T; %[m]
Uo = 1; %[V]
D = 0.019; %Diameter [m]
R = D/2; % Radius [m]
delta_t = .25e-8; %[s]
k = 2*pi*lambda^-1; %N�mero de ondas
K = 1; % Constant of the output voltage

ncycle = 1; %Number of cycles
% ncycle = 17; %Number of cycles

%---% Exercise A %---%
[xvector] = (0.001:0.001:0.200);

%A.1 -- %%%%%%%%%%%%%%%%%%%%%%%%
[yvector] = 0*R;

%A.2 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [yvector] = 0*R;

%---% Exercise B %---%
% [yvector] = (-0.0200:0.0001:0.0200);

%B.1 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.005;

%B.2 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.020;

%B.3 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.060;

%B.4 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.100;

%---% Exercise C %---%
%C.1 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.005;
% [yvector] = 0;

%C.2 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.005;
% [yvector] = 0.003;

%C.3 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.005;
% [yvector] = 0.0095;

%C.4 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.005;
% [yvector] = 0.015;

% C.5 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.020;
% [yvector] = 0;

%C.6 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.020;
% [yvector] = 0.003;

%C.7 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.020;
% [yvector] = 0.0095;

%C.8 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.020;
% [yvector] = 0.015;

%C.9 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.060;
% [yvector] = 0;

%C.10 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.060;
% [yvector] = 0.003;

%C.11 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.060;
% [yvector] = 0.0095;

%C.12 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.060;
% [yvector] = 0.015;

%---% Exercise D %---%
%D.1 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = (0.005:0.0005:0.120);
% [yvector] = (-0.020:0.0001:0.020);

%D.2 -- %%%%%%%%%%%%%%%%%%%%%%%%
%  [xvector] = (0.010:0.0005:0.150);
%  [yvector] = (-0.020:0.0001:0.020);

%---% Exercise E %---%
%E.1 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.005;
% [yvector] = 0;

%E.2 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.005;
% [yvector] = 0.003;

%E.3 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.020;
% [yvector] = 0;

%E.4 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.020;
% [yvector] = 0.003;

%E.5 -- %%%%%%%%%%%%%%%%%%%%%%%%
% [xvector] = 0.060;
% [yvector] = 0;

n = 1;

for N=1:length(yvector)
    for M=1:length(xvector)
        x = xvector(M);
        y = yvector(N);
        [E, P,  Pi, phii, omega, vn, ka, S, t, to, t1, t2, t_vn, t_pconv, t_econv] = CalculatedBeamPressure(x, c, R, y , delta_t, Uo, Fface, rho, T, lambda, ncycle, K);
        [E_n, P_n, Pi_n, phii_n] = Normalization(E, P, Pi, phii, c, delta_t, rho, Uo, K);
        Ppp(N, M) = abs(max(P_n) - min(P_n)); % Peak-to-peak amplitude
        Prms(N, M) = (max(P_n) - min(P_n))/(2*sqrt(2)); % Root mean square amplitude
    end  
end

% Axial pressure amplitude for a baffled circular plane piston.
OnAxialPressure = AxialPressure(k, xvector, R);

%S = 0:S/(length(onAxialPressure)-1):S;

%---% Plotting

%---% Exercise A %---%
% A.1 -- %%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
plot(xvector*1000, OnAxialPressure)
title('Axial pressure amplitude for a baffled circular plane piston')
xlabel('z (mm)')
ylabel('P/(2*rho*c*Uo)')
grid on
grid minor
set(gca,'FontSize', 16);

% figure(4)
% plot(xvector*1000, Prms)
% title('Axial normalized pressure root mean square amplitude (1 cycle)')
% xlabel('x (mm)')
% ylabel('Prms = Ppp/(2*sqrt(2))')
% grid on
% grid minor
% set(gca,'FontSize', 16);

% A.2 -- %%%%%%%%%%%%%%%%%%%%%%%%
% figure(3)
% plot(xvector, OnAxialPressure)
% title('Axial pressure amplitude for a baffled circular plane piston')
% xlabel('x (m)')
% ylabel('P/(2*rho*c*Uo)')

% figure(4)
% plot(xvector, Prms)
% title('Axial normalized pressure root mean square amplitude (30 cycles)')
% xlabel('x (m)')
% ylabel('Prms = Ppp/(2*sqrt(2))')

%---% Exercise B %---%
% B.1 -- %%%%%%%%%%%%%%%%%%%%%%%%
% figure(5)
% plot(yvector*1000, Prms)
% % title('Calculated transmit-receive mode beam plot at 5 mm (1 cicle)')
% title('Calculated transmit-receive mode beam plot at 5 mm (17 cicles)')
% xlabel('Distance off axis (mm)')
% ylabel('Prms = Ppp/(2*sqrt(2))')
% grid on
% grid minor
% set(gca,'FontSize', 16);

% B.2 -- %%%%%%%%%%%%%%%%%%%%%%%%
% figure(6)
% plot(yvector*1000, Prms)
% % title('Calculated transmit-receive mode beam plot at 20 mm (1 cicle)')
% title('Calculated transmit-receive mode beam plot at 20 mm (17 cicles)')
% xlabel('Distance off axis (mm)')
% ylabel('Prms = Ppp/(2*sqrt(2))')
% grid on
% grid minor
% set(gca,'FontSize', 16);

% B.3 -- %%%%%%%%%%%%%%%%%%%%%%%%
% figure(7)
% plot(yvector*1000, Prms)
% % title('Calculated transmit-receive mode beam plot at 60 mm (1 cicle)')
% title('Calculated transmit-receive mode beam plot at 60 mm (17 cicles)')
% xlabel('Distance off axis (mm)')
% ylabel('Prms = Ppp/(2*sqrt(2))')
% grid on
% grid minor
% set(gca,'FontSize', 16);

% B.4 -- %%%%%%%%%%%%%%%%%%%%%%%%
% figure(8)
% plot(yvector*1000, Prms)
% title('Calculated transmit-receive mode beam plot at 100 mm (1 cicle)')
% % title('Calculated transmit-receive mode beam plot at 100 mm (17 cicles)')
% xlabel('Distance off axis (mm)')
% ylabel('Prms = Ppp/(2*sqrt(2))')
% grid on
% grid minor
% set(gca,'FontSize', 16);

%---% Exercise C %---%
% From C.1 to C.7-- %%%%%%%%%%%%%%%%%%%%%%%%

% figure(1)
% subplot(2,2,1)
% plot(t, phii_n, 'r')
% ylabel('phii(r,t)')
% xlabel('Time (s)')
% % title('Normalized velocity potential IR (1 cycle)')
% title('Normalized velocity potential IR (17 cycles)')
% grid on
% grid minor
% set(gca,'FontSize', 14);
% 
% subplot(2,2,2)
% plot(t, Pi_n, 'b')
% ylabel('Pi(r,t)')
% xlabel('Time (s)')
% % title('Normalized pressure IR (1 cycle)')
% title('Normalized pressure IR (17 cycles)')
% grid on
% grid minor
% set(gca,'FontSize', 14);
% 
% 
% subplot(2,2,3)
% plot(t_pconv, P_n, 'm')
% ylabel('P(r,t) = vn(t)*Pi(r,t)')
% xlabel('Time (s)')
% % title('Normalized pressure (1 cycle)')
% title('Normalized pressure (17 cycles)')
% grid on
% grid minor
% set(gca,'FontSize', 14);
% 
% 
% subplot(2,2,4)
% plot(t_econv, E_n, 'k')
% ylabel('E(t)')
% xlabel('Time (s)')
% % title('Normalized output voltage (1 cycle)')
% title('Normalized output voltage (17 cycles)')
% grid on
% grid minor
% set(gca,'FontSize', 14);


%---% Exercise D %---%
%D.1 -- %%%%%%%%%%%%%%%%%%%%%%%%
%  figure(2)
%  pcolor(xvector*1000, yvector*1000, Ppp)
%  xlabel('z(mm)')
%  ylabel('x(mm)')
% %  title('Pressure fields - Maximum peak-to-peak pressure (1 cycle)')
%  title('Pressure fields - Maximum peak-to-peak pressure (17 cycles)')
%  colorbar
%  shading interp
%  grid on
%  grid minor
%  set(gca,'FontSize', 14);

%  figure(9)
%  mesh(xvector*1000, yvector*1000, Ppp)
%  xlabel('z(mm)')
%  ylabel('x(mm)')
%  zlabel('|P(r,t)|')
% % title('Pressure fields - Maximum peak-to-peak pressure (1 cycle)')
%  title('Pressure fields - Maximum peak-to-peak pressure (17 cycles)')
%  colorbar
%  shading interp
%  grid on
%  grid minor
%  set(gca,'FontSize', 14);

%---% Exercise E %---%
% E.1 -- %%%%%%%%%%%%%%%%%%%%%%%%
% figure(10)
% plot(t_econv, E_n, 'k')
% ylabel('E(t)')
% xlabel('Time (s)')
% % title('Normalized output voltage (x = 0 mm and z = 5 mm and 1 cycle)')
% % title('Normalized output voltage (x = 3 mm and z = 5 mm and 1 cycle)')
% % title('Normalized output voltage (x = 0 mm and z = 20 mm and 1 cycle)')
% % title('Normalized output voltage (x = 3 mm and z = 20 mm and 1 cycle)')
% title('Normalized output voltage (x = 0 mm and z = 60 mm and 1 cycle)')
% grid on
% grid minor
% set(gca,'FontSize', 14);