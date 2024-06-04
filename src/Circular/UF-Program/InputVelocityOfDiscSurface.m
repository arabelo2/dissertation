function vn = InputVelocityOfDiscSurface (Uo, f0, t, phase)
    vn = Uo*sin(2*pi*f0*t + phase);       
end
