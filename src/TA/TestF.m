vn=resample(hdf5read('ref_pulse-40MHz.h5', 'ascan'), 100, 40);
% k = numel(vn) 
% vn = vn(find(vn~=0, 1, 'first'):find(vn~=0, 1, 'last'));
% vn(end:end+k - numel(vn))=0;
% size(vn)
% plot(vn)
%vn = vn(1:640);
N = length(vn); % the number of data points
fs = 100e6; % sample rate, number of samples per second
fc = 5e6; 
dt = 1/fs; % time increment (in seconds per sample)
t=dt*(0:N-1); % time domain (in seconds)
figure();
plot(vn);
 
 %spectrum
 VN=fftshift(fft(vn)); % DFT and shift center to zero
 df=fs/N; % the frequency increment
 f=-fs/2:df:fs/2-df; % create the frequency axis
 figure()
 plot(f,abs(VN)/max(abs(VN)), 'ro')
 hold on; plot(f,abs(VN)/max(abs(VN)))
 grid on
 grid minor
 
 % Max frequency
 f(VN == max(VN))
 
hw = hann(find(vn == min(vn)) + 20);
hw(end:end+numel(vn) - numel(hw)) = 0;
plot(hw)
hold on
plot(vn)
 

[bc,ac] = butter(8,fc/(fs/2));
vn_filtered = filter(bc,ac,vn);
hold on
hold on; plot(vn_filtered)

 %spectrum
 VN=fftshift(fft(vn_filtered)); % DFT and shift center to zero
 df=fs/N; % the frequency increment
 f=-fs/2:df:fs/2-df; % create the frequency axis
 hold on; plot(f,abs(VN)/max(abs(VN)), 'ro')
 hold on; plot(f,abs(VN)/max(abs(VN)))
 grid on
 grid minor