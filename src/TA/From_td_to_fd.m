N = length(vn); % the number of data points
fs = 100e6; % sample rate, number of samples per second
dt = 1/fs; % time increment (in seconds per sample)
t=dt*(0:N-1); % time domain (in seconds)
fc = 5e6;  % carrier frequency
%x=sin(2*pi*fc*t); % sine
%vn = x;
plot(vn);
figure();
plot(t,vn);
plot(t*1500,vn);

%spectrum
VN=fftshift(fft(vn)); % DFT and shift center to zero
df=fs/N; % the frequency increment
f=-fs/2:df:fs/2-df; % create the frequency axis
figure()
plot(f,abs(VN))
grid on
grid minor

% Max frequency
f(VN == max(VN))

% reverse procedure
% y=ifft(VN); % inverse fft
% N2=length(VN); % determine the length of the signal
% dt2=1/fs; % determine the time increment
% tim=0:dt2:(N2-1)*dt2; %create the time axis
% figure()
% plot(tim,y)

%
% https://www.mathworks.com/matlabcentral/answers/38552-from-time-domain-to-frequency-domain-and-back-to-time-domain
%