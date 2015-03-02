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
        if ~( Eseq((i-1)*24+j) >= 0 && Eseq((i-1)*24+j) < 2000)
            Eseq((i-1)*24+j) = 0;
        end
    end
end
figure;
plot(Eseq,'.-')
A = zeros(365,1);
B = zeros(365,1);
Year27 = Eseq;
figure;
plot(Year27,'.-')
Year = reshape(Year27,365*24,27);
Day = reshape(Year27,24,365*27);
Day = sum(Day);
Day = Day';
figure;
plot(Day)

Hournewnozero = zeros(size(Year27));
count = 0;
for i = 1:365*27*24
    if Year27(i) > 10
        count = count + 1;
        Hournewnozero(count) = Year27(i);        
    end
end
Hournewnozero = Hournewnozero(1:count);
[h,center] = hist(Hournewnozero,100);
h = h/sum(h);
figure;
bar(center,h);
xlabel('Harvested Energy in One Hour');
ylabel('PDF');
grid on;
acch = zeros(size(h));
acch(1) = h(1);
for i = 2:length(h)
    acch(i) = acch(i-1) + h(i);
end
figure;
plot(center,acch,'.-');
xlabel('Harvested Energy in One Hour');
ylabel('CDF');
grid on;

Yearmax = max(Year,[],2);
figure;
plot(Yearmax,'.-')
xlabel('Hour');
ylabel('Maximal harvested energy');
% for i = 1+1:365*24-1
%     if Yearmax(i) > Yearmax(i-1)+500 && Yearmax(i) > 1.5 *Yearmax(i-1)
%         Yearmax(i) = (Yearmax(i-1) + Yearmax(i+1))/2;
%     end
% end
% syms a b
% for i = 1:365
%     Day = Yearmax((i-1)*24+1:i*24);
%     m = max(Day);
%     s = sum(Day);
%     fz = @(a) 2*(a*12/pi*sqrt(1-((m-a)/a)^2)+ (m-a)*12*acos(-(m-a)/a)) - s;
%     A(i) = fzero(fz,m);
%     B(i) = m - A(i);
% end
% A
% B
% figure;
% plot(A)
% figure;
% plot(B)
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
Sind = zeros(365,1);
for d = 1:365
    f = @(x) E0(d)*max(sin(lat)* sin(Theta(d)) + cos(lat)*cos(Theta(d))*cos(pi*(x-12)/12),0);
    for h = 1:24
        Sinh((d-1)*24+h) = integral(f,h-1 + LAT(d),h + LAT(d));        
    end
    Sind(d) = integral(f,0+ LAT(d),24+ LAT(d));
end
figure;
plot(Sinh,'.-');
xlabel('Hour');
ylabel('Energy that can be harvested');

L = length(Sinh);
NFFT = 2^nextpow2(L);
y = fft(Sinh,NFFT)/L;
Fs = 1;
f = Fs/2*linspace(0,1,NFFT/2+1);
figure;
plot(f,2*abs(y(1:NFFT/2+1)),'g') 
title('Single-Sided Amplitude Spectrum of y3(t)')
xlabel('Frequency (1/h)')
ylabel('|Y3(f)|')
grid on

Yearnew = zeros(size(Year));
Yearnew1 = zeros(size(Year));
Yearnew2 = zeros(size(Year));
NorSinh = zeros(size(Sinh));
for d = 1:365
    xx = Sinh((d-1)*24+1:d*24);
    NorSinh((d-1)*24+1:d*24) = xx/norm(xx);
end
Daynew = zeros(size(Day));
for i = 1:27
    Yearnew(:,i) = Year(:,i)./Sinh;
    Yearnew2(:,i) = Year(:,i)./Yearmax;    
    for d = 1:365
        Co = sum(Year((d-1)*24+1:d*24,i).*NorSinh((d-1)*24+1:d*24));
        Yearnew1((d-1)*24+1:d*24,i) = Year((d-1)*24+1:d*24,i) - Co*NorSinh((d-1)*24+1:d*24);
    end
    Daynew((i-1)*365+1:i*365) = Day((i-1)*365+1:i*365)./Sind;
end
figure;
plot(Daynew,'.-r');
xlabel('Day');
ylabel('E(d)/E_{max}(d)');

Daynewnozero = zeros(size(Daynew));
count = 0;
for i = 1:365*27
    if Daynew(i) > 10
        count = count + 1;
        Daynewnozero(count) = Daynew(i);        
    end
end
Daynewnozero = Daynewnozero(1:count);
[h,center] = hist(Daynewnozero,100);
h = h/sum(h);
figure;
bar(center,h);
xlabel('Harvested Energy in One Day');
ylabel('PDF');
grid on;

acch = zeros(size(h));
acch(1) = h(1);
for i = 2:length(h)
    acch(i) = acch(i-1) + h(i);
end
figure;
plot(center,acch,'.-');
xlabel('Harvested Energy in One Day');
ylabel('CDF');
grid on;

Year27new = reshape(Yearnew,27*365*24,1);
Year27new1 = reshape(Yearnew1,27*365*24,1);
Year27new2 = reshape(Yearnew2,27*365*24,1);
for i = 1:27*365*24
    if  Year27new(i) > 2e3
        Year27new(i) = 0;
    end
    if ~(Year27new2(i) >= 0 && Year27new2(i) <= 1)
        Year27new2(i) = 0;
    end
end

figure;
plot(Year27new,'.-')
xlabel('Hour')
ylabel('E(t)/E_{max}(t)')

figure;
plot(Year27new2,'.b')
xlabel('Hour')
ylabel('E(t)/E_{max1}(t)')

figure;
plot(Year27new1,'.g')
xlabel('Hour')
ylabel('E(t) - E(t)(E(t)* v(t))')

L = length(Year27new);
NFFT = 2^nextpow2(L);
y = fft(Year27new,NFFT)/L;
Fs = 1;
f = Fs/2*linspace(0,1,NFFT/2+1);
figure;
plot(f,2*abs(y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y1(t)')
xlabel('Frequency (1/h)')
ylabel('|Y1(f)|')
grid on

y = fft(Year27new1,NFFT)/L;
Fs = 1;
f = Fs/2*linspace(0,1,NFFT/2+1);
figure;
plot(f,2*abs(y(1:NFFT/2+1)),'g') 
title('Single-Sided Amplitude Spectrum of y2(t)')
xlabel('Frequency (1/h)')
ylabel('|Y2(f)|')
grid on


