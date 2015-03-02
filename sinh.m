function y = sinh(x)
    n = 1:365;
    ThetaC = 23.45* sin(2*pi*(284+n)/365)*pi/180;
    y = max (sin(lat)* sin(ThetaC(d)) + cos(lat)*cos(ThetaC(d))*cos(pi*(x-12)/12),0);
end