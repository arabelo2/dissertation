N=9000; %numero de pontos
py=120;
s2=load(['E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Temp\PPGEM\Dissertação\Programa\Matuda\0graus-2015-03-01\0graus\DADOS0H0.DAT']); % ou 78, ou 76, veja se tem muita diferença
s2_temp=reshape(s2,N,py);
vn = resample(s2_temp(:, 60), 100, 100);
vn = vn(find(vn == max(vn))-25:find(vn == max(vn)) + 215)/max(abs(vn)); %100MHz
% vn = vn(find(vn == max(vn))-80:find(vn == max(vn)) + 333)/max(abs(vn)); %320MHz
% vn = vn(find(vn == max(vn)):find(vn == max(vn)) + 2*1024)/max(abs(vn)); %320MHz
save sinalvn0h0 vn
figure(1)
subplot(2,1,1)
plot(vn)
grid on
grid minor
hold on
legend('DADOS0H0')
% From time domain to frequency domain
N = length(vn); % the number of data points
fs = 100e6; % sample rate, number of samples per second
dt = 1/fs; % time increment (in seconds per sample)
t=dt*(0:N-1); % time domain (in seconds)
fc = 5e6;  % carrier frequency
%spectrum
VN=fftshift(fft(vn)); % DFT and shift center to zero
df=fs/N; % the frequency increment
f=-fs/2:df:fs/2-df; % create the frequency axis
figure(1);
subplot(2,1,2)
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