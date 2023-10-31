function vn = InputVelocityOfDiscSurface (Uo, f0, t_vn)
    vn = Uo*sin(2*pi*f0*t_vn);       
end

% while (2*pi*f0*t(m) <= D)
%     vn(m) = Uo*exp(j*2*pi*f0*t(m)); % vn(m) = Uo*(cos(2*pi*f0*t(m)) + j*sen(2*pi*f0*t(m)));  
%     m = m + 1;
% end
