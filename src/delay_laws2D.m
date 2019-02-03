function td = delay_laws2D (M, s, Phi, F, c)
Mb = (M -1)/2;
m = 1:M;
em = s*((m - 1) - Mb); % Location of centroids of elements
switch(F)
    % Steering only case
    case inf
        if Phi > 0
            td = s*sind(Phi)*(m - 1)/c;
        else
            td = s*sind(abs(Phi))*(M - m)/c;
        end
        % Steering and focusing case
    otherwise 
        r1 = sqrt(F^2 + (Mb*s)^2 + 2*F*Mb*s*sind(Phi));
        rm = sqrt(F^2 + em.^2 - 2*F*em*sind(Phi));
        rM = sqrt(F^2 + (Mb*s)^2 + 2*F*Mb*s*sind(abs(Phi)));
        if Phi > 0
            td = (r1 - rm)/c;
        else
            td = (rM - rm)/c;
        end
end
    
% Reference
% Jr, Lester W. Schmerr. Fundamentals of Ultrasonic Phased Arrays. Solid Mechanics and Its Applications. Springer International Publishing, 2015. //www.springer.com/us/book/9783319072715.