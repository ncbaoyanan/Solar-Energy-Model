clear all
% close all
clc
lat = 39.742*pi/180;
lon = 105.18;
D = 1:365;
Dangle = 2*pi*(D-1)/365;
E0 = 1.000110 + 0.034221*cos(Dangle) + 0.001280*sin(Dangle) + 0.000719*cos(2*Dangle) + 0.000077*sin(2*Dangle);

Theta = 0.006918 - 0.399912*cos(Dangle) + 0.070257*sin(Dangle) - 0.006758*cos(2*Dangle) + 0.000907*sin(2*Dangle) - 0.002697*cos(3*Dangle) + 0.00148*sin(3*Dangle);

Et = (0.000075 + 0.001868*cos(Dangle) - 0.032077*sin(Dangle) - 0.014615*cos(2*Dangle) - 0.04089*sin(2*Dangle))*229.18/60;

LAT = 4*(105-lon)/60 + Et;

d=[31 28 31 30 31 30 31 31 30 31 30 31];
beginday=[1 32 60    91   121   152   182   213   244   274   305   335];
endday=[31 59 90 120 151 181 212 243 273 304 334 365];

height = 1.829;
a0 = 0.4237 - 0.00821*(6-height)^2;
a1 = 0.5055 + 0.00595*(6.5-height)^2;
k = 0.2711 + 0.01858*(2.5-height)^2;

% n = 1:365;
%     ThetaC = 23.45* sin(2*pi*(284+n)/365)*pi/180;

A = [1202 1187 1164 1130 1106 1092 1093 1107 1136 1136 1190 1204];
B = [0.141 0.142 0.149 0.164 0.177 0.185 0.186 0.182 0.165 0.152 0.144 0.141];
C = [0.103 0.104 0.109 0.120 0.130 0.137 0.138 0.134 0.121 0.111 0.106 0.103];

Eh = zeros(365*24,1);
m = 1;
for d = 1:365
    
%     if d > endday(m)
%         m = m + 1;
%     end
%     f = @(x) E0(d)*A(m)*C(m)*(exp(-B(m)./(max(sin(lat)* sin(Theta(d)) + cos(lat)*cos(Theta(d))*cos(pi*(x-12)/12),0)))).*max(sin(lat)* sin(Theta(d)) + cos(lat)*cos(Theta(d))*cos(pi*(x-12)/12),0);
%     
%     cosh = - tan(lat)*tan(Theta(d));
%     hss  = acos(cosh)*180/pi;
%     alpha = 0.409 + 0.5016*sin(hss- 60);
%     beta = 0.6609 -0.4767*sin(hss - 60);
    
%     f = @(x) E0(d).*max(sin(lat)* sin(Theta(d)) + cos(lat)*cos(Theta(d))*cos(pi*(x-12)/12),0);
    f = @(x) E0(d)*(a0+a1*exp(-k./(max(sin(lat)* sin(Theta(d)) + cos(lat)*cos(Theta(d))*cos(pi*(x-12)/12),0)))).*max(sin(lat)* sin(Theta(d)) + cos(lat)*cos(Theta(d))*cos(pi*(x-12)/12),0);
    for h = 1:24
        Eh((d-1)*24 + h) = integral(f,h-1 + LAT(d),h + LAT(d));        
    end  
%     total = integral(f,0+ LAT(d),24+ LAT(d));
%     for h = 1:24    
%     Eh((d-1)*24 + h) = total*pi/24*(alpha+ beta*cos(pi*(h-0.5+ LAT(d)-12)/12)).*(cos(pi*(h-0.5+ LAT(d)-12)/12) - cos(hss))./(sin(hss) - (2*pi*hss/360)*cos(hss));
%     end
end
figure;
plot(Eh,'.-');
xlabel('Hour');
ylabel('Energy that can be harvested');