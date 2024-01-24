
%% Video of a moving pressure

% Parameters
close all;
run_c;

% Create the writer object and open it
writerObj = VideoWriter('AxialPressureVideo x = 0 mm and 1 cycles', 'MPEG-4'); % Name of video file
open(writerObj);

% Create plotting figure
fig1 = figure(1);

% Loop over all the points

for indx = 1:length(x)
    for indz = 1:length(z)
        legendInfo = ['x = ',num2str(x(indx)), ' | z = ', num2str(z(indz))];
        plot(t_conv_c{indx, indz}*c1*1e3, P_c{indx, indz}) % P_c{yy, xx} = rho*conv(h_temp, diff(vn)/(t_temp(2) - t_temp(1)));
        title(['Axial Pressure amplitude for a baffled circular plane piston with ', num2str(nc),' cycle(s)'], 'Color', 'k', 'interpreter', 'latex')
        xlabel('c.t(mm)')
        ylabel('Pressure')
        legend(legendInfo)
        axis([0 200 -11e13 11e13])
        grid on
        grid minor
        pause(0.05)
        
        % Retrieve the current frame of the figure
        F = getframe(fig1);
        
        % Write the current frame to the writer object
        writeVideo(writerObj, F);
    end
end

% Close the writer object. File appears in current folder
close(writerObj);
disp('Video file written successfully!')