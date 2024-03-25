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
lambda = c1 / f0; % [m]
STEP = lambda / 10;
NF = R^2/lambda;% Near Field Length or Transition from Near Field to Far Field
k = 2*pi/lambda; % wave number
Uo = 1; % [V]
K = 1; % Constant of the output voltage
nc = 27; % Number of cycles
Dhydrophone = 0.6e-3;
discretization = 51;

xmin = 0.002;  % z-axis
xmax = 0.102; % z-axis

ymin = 0; % x-axis
ymax = 0; % x-axis

xpoints = 667;
ypoints = 1;

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
        [t_temp, h_temp] = vpirOfCircularPistonlikeTransducersWScanner(x(xx), y(yy), D, c1, fs);
        h_c{yy, xx} = h_temp;
        t_c{yy, xx} = t_temp;
        
        %%%% texcitation{xx, yy, zz} = texcitation_temp;
        v_c{yy, xx} = vn;

        % Pressure distribution
        Pi_temp = rho*diff(h_temp)/(t_temp(2) - t_temp(1)); 
        Pi_c{yy, xx} = Pi_temp;
            
        % Pressure waveform
        % Transient pressure
        p_temp = rho*conv(h_temp, diff(vn));
        
        P_c{yy, xx} = p_temp;
        t_conv_temp = t_temp(1) + (t_temp(2) - t_temp(1))*(0:1:length(p_temp)-1);
        t_conv_c{yy, xx} = t_conv_temp;

        N = round(fs/f0); % Sample
        mid = floor(length(p_temp)/2);
        index = mid - N : 1 : mid + N; % garantindo 2 ciclos para buscar max(abs(p))
        p_temp = p_temp(index); % vetor das pressões somente da parte central (2 ciclos) da onda completa
            
        % Peak amplitude
        Pp_c(yy, xx) = max(abs(p_temp));
    end
end

% Axial pressure amplitude for a baffled circular plane piston.
OnAxialPressure = AxialPressure(k, x, R);

z = x;
x = y;

Pp_c_prg = Pp_c;

profile viewer 
toc;
