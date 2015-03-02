% clear all
close all
clc
load('dat.mat');
value=dat;

[lx ly] = size(value);

Validh = ones(365*24*27,1);
Eseq = zeros(lx*ly,1);
for i = 1:lx
    for j = 1:ly
        Eseq((i-1)*24+j) = value(i,j);
        if ~( Eseq((i-1)*24+j) >= 0 && Eseq((i-1)*24+j) < 2000)
            Eseq((i-1)*24+j) = 0;
            Validh((i-1)*24+j) = 0;
        end
    end
end
hour = Eseq;

Isc = 1367;         % 假设最大辐射强度不超过Isc 单位是瓦特/平方米
lat = 39.742*pi/180;    % 纬度
lon = 105.18;           % 经度
D = 1:365;              % 天
Dangle = 2*pi*(D-1)/365;  % 跟地球运转有关的角度
% E0是日地距离
E0 = 1.000110 + 0.034221*cos(Dangle) + 0.001280*sin(Dangle) + 0.000719*cos(2*Dangle) + 0.000077*sin(2*Dangle);
% Theta 
Theta = 0.006918 - 0.399912*cos(Dangle) + 0.070257*sin(Dangle) - 0.006758*cos(2*Dangle) + 0.000907*sin(2*Dangle) - 0.002697*cos(3*Dangle) + 0.00148*sin(3*Dangle);
% 一天长度的变化，和纬度有关
Et = (0.000075 + 0.001868*cos(Dangle) - 0.032077*sin(Dangle) - 0.014615*cos(2*Dangle) - 0.04089*sin(2*Dangle))*229.18/60;
LAT = 4*(105-lon)/60 + Et;
height = 1.829;
a0 = 0.4237 - 0.00821*(6-height)^2;
a1 = 0.5055 + 0.00595*(6.5-height)^2;
k = 0.2711 + 0.01858*(2.5-height)^2;

Eh = zeros(365*24,1);
m = 1;
for d = 1:365
%       f = @(x) Isc*E0(d)*(a0+a1*exp(-k./(max(sin(lat)* sin(Theta(d)) + cos(lat)*cos(Theta(d))*cos(pi*(x-12)/12),0)))).*max(sin(lat)* sin(Theta(d)) + cos(lat)*cos(Theta(d))*cos(pi*(x-12)/12),0);
    f = @(x) Isc*E0(d)*(0.2710 +(1-0.2939)*(a0+a1*exp(-k./(max(sin(lat)* sin(Theta(d)) + cos(lat)*cos(Theta(d))*cos(pi*(x-12)/12),0))))).*max(sin(lat)* sin(Theta(d)) + cos(lat)*cos(Theta(d))*cos(pi*(x-12)/12),0);
%     f = @(x) Isc*E0(d).*max(sin(lat)* sin(Theta(d)) + cos(lat)*cos(Theta(d))*cos(pi*(x-12)/12),0);
    for h = 1:24
        Eh((d-1)*24 + h) = integral(f,h-1 + LAT(d),h + LAT(d));        
    end
end

Len = 365*24*27;
Datagood = zeros(Len,1);
Datanozero = zeros(Len,1);
count = 0;
for i = 1:Len
    if Eh(mod(i-1,365*24)+1) > Isc*0.05
        Datagood(i) = hour(i)/Eh(mod(i-1,365*24)+1);
             
%         if Datagood(i) > 1
%              Datagood(i) = 1;
%         end
        count = count + 1;   
        Datanozero(count) = Datagood(i);
    else
        Datagood(i) = 0;
    end
end

istart = (26)*365*24 + 1;
iend = istart + 31*24 - 1;

d201001 = Datagood(istart:1:iend);
figure;
plot(d201001,'.-');
figure;
plot(Validh(istart:1:iend),'.-');
save d201001.mat d201001