function OnAxialPressure = AxialPressure(k, xvector, R)
OnAxialPressure = abs(sin((1/2).*k*xvector.*(sqrt(1 + R^2./xvector.^2) -1)));