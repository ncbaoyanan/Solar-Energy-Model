clear all
close all
clc
%dat=[x0911(:,2:2:48); x9508_star(:,2:2:48)];
Cos = zeros(24,1);
f = @(x) cos(pi*(x-12)/12);
for i = 1:24
    Cos(i) = quad( f,i-1,i);
end
save Cos
load Cos
load('dat.mat');
value=dat;

[lx ly] = size(value);

Eseq = zeros(lx*ly,1);
for i = 1:lx
    for j = 1:ly
        Eseq((i-1)*24+j) = value(i,j);
        if ~( Eseq((i-1)*24+j) >= 0.1 && Eseq((i-1)*24+j) < 2000)
            Eseq((i-1)*24+j) = 0.1;
        end
    end
end
Year27log = log(Eseq);
figure;
plot(Year27log,'.-')

L = length(Year27log);
NFFT = 2^nextpow2(L);
y = fft(Year27log,NFFT)/L;
Fs = 1;
f = Fs/2*linspace(0,1,NFFT/2+1);
figure;
plot(f,2*abs(y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y1(t)')
xlabel('Frequency (1/h)')
ylabel('|Y1(f)|')
grid on


lat = 39.742*pi/180;
lon = 105.18;
D = 1:365;
Dangle = 2*pi*(D-1)/365;
E0 = 1.000110 + 0.034221*cos(Dangle) + 0.001280*sin(Dangle) + 0.000719*cos(2*Dangle) + 0.000077*sin(2*Dangle);
Theta = 0.006918 - 0.399912*cos(Dangle) + 0.070257*sin(Dangle) - 0.006758*cos(2*Dangle) + 0.000907*sin(2*Dangle) - 0.002697*cos(3*Dangle) + 0.00148*sin(3*Dangle);
Et = (0.000075 + 0.001868*cos(Dangle) - 0.032077*sin(Dangle) - 0.014615*cos(2*Dangle) - 0.04089*sin(2*Dangle))*229.18/60;
LAT = 4*(105-lon)/60 + Et;

n = 1:365;
ThetaC = 23.45* sin(2*pi*(284+n)/365)*pi/180;
Alpha = pi*(n-12)/12;
Sinh = zeros(365*24,1);
Sinhh = zeros(365*24,1);
Sind = zeros(365,1);
for d = 1:365
    f = @(x) 1367*E0(d)*max(sin(lat)* sin(Theta(d)) + cos(lat)*cos(Theta(d))*cos(pi*(x-12)/12),0);
    for h = 1:24
        Sinh((d-1)*24+h) = integral(f,h-1 + LAT(d),h + LAT(d));   
        if Sinh((d-1)*24+h) < 0.1
            Sinh((d-1)*24+h) = 0.1;
        else
            Sinhh((d-1)*24+h) = Sinh((d-1)*24+h)/1367/E0(d);
        end       
    end    
    Sind(d) = integral(f,0+ LAT(d),24+ LAT(d));
end

Sinh = log(Sinh);
Year27log = reshape(Year27log,365*24,27);

Pure = zeros(size(Year27log));
for i = 1:27
    Pure(:,i) = (Year27log(:,i) - Sinh).*Sinhh;
end
Pure = reshape(Pure,365*24*27,1);
figure;
plot(Pure,'.-r')


% Yearnew1 = zeros(size(Year27log));
% NorSinh = zeros(size(Sinh));
% for d = 1:365
%     xx = Sinh((d-1)*24+1:d*24);
%     xx = xx - mean(xx);
%     NorSinh((d-1)*24+1:d*24) = xx/norm(xx);
% end
% figure;
% plot(NorSinh,'.-')
% for i = 1:27
%     for d = 1:365
%         yy= Year27log((d-1)*24+1:d*24,i) - mean(Year27log((d-1)*24+1:d*24,i));
%         Co = sum(yy.*NorSinh((d-1)*24+1:d*24));
%         Yearnew1((d-1)*24+1:d*24,i) = yy - Co*NorSinh((d-1)*24+1:d*24);
%     end
% end
% Year27new1 = reshape(Yearnew1,27*365*24,1);
% figure;
% plot(Year27new1,'.-g')
% xlabel('Hour')
% ylabel('E(t) - E(t)(E(t)* v(t))')
% 
% L = length(Year27new1);
% NFFT = 2^nextpow2(L);
% y = fft(Year27new1,NFFT)/L;
% Fs = 1;
% f = Fs/2*linspace(0,1,NFFT/2+1);
% figure;
% plot(f,2*abs(y(1:NFFT/2+1))) 
% title('Single-Sided Amplitude Spectrum of y2(t)')
% xlabel('Frequency (1/h)')
% ylabel('|Y2(f)|')
% grid on
% 
% 
