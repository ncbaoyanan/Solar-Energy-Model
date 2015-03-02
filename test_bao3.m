clear all
close all
clc
height = 1.829;
a0 = 0.4237 - 0.00821*(6-height)^2;
a1 = 0.5055 + 0.00595*(6.5-height)^2;
k = 0.2711 + 0.01858*(2.5-height)^2;
lat = 39.742*pi/180;
n = 1:365;
    ThetaC = 23.45* sin(2*pi*(284+n)/365)*pi/180;
n = 1:24;
Eh = zeros(365*24,1);
for d = 1:365
    f = @(x) (1 + 0.034*cos(2*pi*d/365))*(a0+a1*exp(-k./(max(sin(lat)* sin(ThetaC(d)) + cos(lat)*cos(ThetaC(d))*cos(pi*(x-12)/12),0)))).*max(sin(lat)* sin(ThetaC(d)) + cos(lat)*cos(ThetaC(d))*cos(pi*(x-12)/12),0);
    for h = 1:24
        Eh((d-1)*24+h) = integral(f,h-1,h);
    end
end
figure;
plot(Eh,'.');
xlabel('Hour');
ylabel('Energy that can be harvested');

load('dat.mat');
value=dat;
Ehh = reshape(Eh,24,365);
Ehh = Ehh';
figure;
plot(Ehh);

[lx ly] = size(value);

Eseq = zeros(lx*ly,1);
for i = 1:lx
    for j = 1:ly
        Eseq((i-1)*24+j) = value(i,j);
        if ~( Eseq((i-1)*24+j) >= 0 && Eseq((i-1)*24+j) < 2000)
            Eseq((i-1)*24+j) = 0;
        end
    end
end

Year27 = Eseq;
figure;
plot(Year27,'.')
xlabel('Hour')
ylabel('E(t)')
Year = reshape(Year27,365*24,27);
Yearnew = zeros(size(Year));
for i = 1:27
    Yearnew(:,i) = Year(:,i)./Eh;
end

Year27new = reshape(Yearnew,27*365*24,1);

for i = 1:27*365*24
    if  ~(Year27new(i)>=0 && Year27new(i)<= 2000)
        Year27new(i) = 0;
    end
end

figure;
plot(Year27new,'.')
xlabel('Hour')
ylabel('E(t)/E_{max}(t)')