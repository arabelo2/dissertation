close all

% Adds the specified folders to the top of the search path for the current MATLAB® session.
addpath('E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\code\src\Array')
addpath('E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\code\src\Rectangular')
addpath('E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Temp\PPGEM\Dissertação\Programa\Matuda')

% analyze the frequency spectrum of a time-based data
fs_original = 1.920030720491530e+07;
fs_original = ceil(fs_original);

% New desired sampling frequency
fs = 100*fs_original;

% Calculate the greatest common divisor (GCD) for p and q
gcdVal = gcd(fs_original, fs);

% Calculate p and q based on the GCD
p = fs / gcdVal;
q = fs_original / gcdVal;

% Load your data: Import your time-domain data into MATLAB.
% v_temp = hdf5read('ref_pulse-40MHz.h5', 'ascan');
% v_temp = resample(v_temp, Fs*1e-6, 40);
for xx = 1:length(z)
    for yy = 1:length(x)
data_i = V{yy, xx};
data_o = P_c{yy, xx};

% Resample the data
data_i_resampled = resample(data_i, p, q);
% data_i_resampled = data_o;

% Apply the FFT: Use the fft function to compute the Fourier transform of 
% your data.
Y = fft(data_i_resampled);
H = fft([h_c{yy, xx}; zeros(length(data_i_resampled) - length(h_c{yy, xx}), 1)]);
X_estimated = ifft(Y ./ H);

% figure(1)
% plot(X_estimated)

         cutoff_freq = 0.4e0+6; % Cutoff frequency (Hz)
         filter_order = 6; % Filter order
         lp_filter = designfilt('lowpassfir', 'FilterOrder', filter_order, 'CutoffFrequency', cutoff_freq, 'SampleRate', fs);

         % Apply the filter to the noisy signal
         filtered_signal = filter(lp_filter, X_estimated);
         filtered_signal = detrend(filtered_signal);
         v_c{yy, xx} = filtered_signal;

                  % Plot the original and filtered signals
         figure(70);
         t_orig = 1/fs*(0:length(X_estimated) - 1);
         plot(t_orig, X_estimated/max(X_estimated), 'b', 'LineWidth', 2);
         ylim([-1 1])
         hold on;
         plot(t_orig, filtered_signal/max(filtered_signal), 'r', 'LineWidth', .5);
         xlabel('Time');
         ylabel('Signal');
         title('Estimated Signal vs. Filtered Signal');
         legend('Estimated Signal', 'Filtered Signal');
         grid on;
         grid minor;
         hold off;  

% Compute the two-sided spectrum P2 and the single-sided spectrum P1
P2 = abs(Y/length(data_i_resampled));
P1 = P2(1:floor(length(data_i_resampled)/2)+1);
P1(2:end-1) = 2*P1(2:end-1);

% Define the frequency domain f based on the new sampling frequency
f = (0:(length(data_i_resampled)/2))*fs/length(data_i_resampled);

% Plot the single-sided amplitude spectrum
figure(2)
plot(f, P1)
%xlim([0 1]*1e+07)
title('Single-Sided Amplitude Spectrum of Resampled Data')
xlabel('Frequency (Hz)')
ylabel('|P1(f)|')

pause(.5)
    end
end

% Here, yourData is the variable containing your time-domain data, and Fs
% is the sampling frequency. The fft function computes the discrete Fourier
% transform of the signal, and the power spectrum is plotted as a function
% of frequency. The plot function is then used to visualize the frequency
% spectrum1.

% Remember to replace yourData and Fs with your actual data and sampling
% rate. This will give you a visual representation of the frequency
% components present in your data set.