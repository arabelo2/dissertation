
%% Video of a moving pressure

% Adds the specified folders to the top of the search path for the current MATLAB® session.
addpath('E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\code\src\Array\')
addpath('E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\code\src\Rectangular\')
addpath('E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\code\src\Circular\UF-Program\')
addpath('E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Users\rabeloal\Documents\PPGEM\PMR5234\Program\code\src\Circular\')
addpath('E:\FileHistory\arabelo@hpe.com\RABELOAL11\Data\C\Temp\PPGEM\Dissertação\Programa\Matuda\')

% Parameters
close all;
run_c;

% Create the writer object and open it
writerObj = VideoWriter('AxialPressureVideo x = 0 mm and 17 cycles', 'MPEG-4'); % Name of video file
%writerObj = VideoWriter('AxialPressureVideo x = 0 mm and 5 cycles', 'Motion JPEG AVI'); % Name of video file
writerObj.Quality = 50;

open(writerObj);

% Create plotting figure
fig1 = figure(1);

% Loop over all the points

for indx = 1:length(x)
    for indz = 1:length(z)
         legendInfo = ['x = ',num2str(x(indx)), ' | z = ', num2str(z(indz))];
%         plot(t_conv_c{indx, indz}*c1*1e3, P_c{indx, indz}) % P_c{yy, xx} = rho*conv(h_temp, diff(vn)/(t_temp(2) - t_temp(1)));
%         title(['Axial Pressure amplitude for a baffled circular plane piston with ', num2str(nc),' cycle(s)'], 'Color', 'k', 'interpreter', 'latex')
%         xlabel('c.t(mm)')
%         ylabel('Pressure')
%         legend(legendInfo)
%         axis([0 125 -12e13 12e13])
%         grid on
%         grid minor
%         pause(1)
        
        % Retrieve the current frame of the figure
        % F = getframe(fig1);
        
        % Write the current frame to the writer object
        % writeVideo(writerObj, F);

%         x2 = P_c{indx, indz};
%         [p1,f1,t1] = pspectrum(x2,fs,'spectrogram');
%         waterfall(f1,t1,p1')
%         xlabel('Frequency (Hz)')
%         ylabel('Time (seconds)')
%         wtf = gca;
%         wtf.XDir = 'reverse';
%         view([45 45])
%         pause(1)

          wgt = P_c{indx, indz};          
          [TF,S1,S2] = ischange(wgt, "linear");
          mid = floor(length(S1)/2);
          findex = findchangepts(diff(S1(1:mid)),'Statistic','rms');
          lindex = mid + findchangepts(diff(S1(mid+1:end)),'Statistic','rms');

          subplot(2,1,1); 
          %figure(1)
          plot(t_conv_c{indx, indz}*c1*1e3, wgt, '-')
%           hold on          
%           plot(t_conv_c{indx, indz}*c1*1e3, wgt,'o')          
          xlabel('c.t(mm)')
          % ylim([-12e13, 12e13])
          axis([0 150 -1e14 1e14])
          title({['Axial Pressure amplitude for a baffled circular plane piston with ', num2str(nc),' cycle(s)'], ['Transient and harmonic waveforms']}, 'Color', 'k', 'interpreter', 'latex')
          ylabel('Pressure')
          legend(legendInfo)
          grid on
          grid minor
          hold off
 
          subplot(2,1,2);
          %figure(2)
          wgt(1:findex) = 0;
          wgt(lindex:end) = 0;          
          plot(t_conv_c{indx, indz}*c1*1e3, wgt, 'r-')
%           hold on
%           plot(t_conv_c{indx, indz}*c1*1e3, wgt,'o')
          xlabel('c.t(mm)')
          % ylim([-12e13, 12e13])
          axis([0 150 -1e14 1e14])
          title({['Axial Pressure amplitude for a baffled circular plane piston with ', num2str(nc),' cycle(s)'], ['Only harmonic waveforms']}, 'Color', 'k', 'interpreter', 'latex')
          ylabel('Pressure')
          legend(legendInfo) 
          grid on
          grid minor
          hold off  

          pause(0.1)

        % Retrieve the current frame of the figure
        F = getframe(fig1);
        
        % Write the current frame to the writer object
        writeVideo(writerObj, F);
    end
end

% Close the writer object. File appears in current folder
close(writerObj);
disp('Video file written successfully!')