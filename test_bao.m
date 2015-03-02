clear all
close all
clc
%dat=[x0911(:,2:2:48); x9508_star(:,2:2:48)];
load('dat.mat');
value=dat;

figure;
plot(value);
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
xlabel('Hour');
ylabel('Harvested Energy');

%Eseq = Eseq(1:1000);
figure;
plot(Eseq,'.-');
L = length(Eseq);
NFFT = 2^nextpow2(L);
y = fft(Eseq,NFFT)/L;
Fs = 1;
f = Fs/2*linspace(0,1,NFFT/2+1);
figure;
plot(f,2*abs(y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (1/h)')
ylabel('|Y(f)|')
grid on

Cycle = zeros(365*24,27);
for i = 1:27
    Cycle(:,i) = Eseq((i-1)*365*24+1:i*365*24);
end
year = median(Cycle,2);
figure;
plot(year)

Eseq1 = Eseq;
for i= 1:27
    Eseq1((i-1)*365*24+1:i*365*24) = Eseq1((i-1)*365*24+1:i*365*24) ./ year;
end

for i = 1:lx
    for j = 1:ly
        if ~( Eseq1((i-1)*24+j) >= 0 && Eseq1((i-1)*24+j) < 2000 )
            Eseq1((i-1)*24+j) = 0;
        end
    end
end
figure;
plot(Eseq1)

y = fft(Eseq1,NFFT)/L;
figure;
plot(f,2*abs(y(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y1(t)')
xlabel('Frequency (1/h)')
ylabel('|Y1(f)|')
grid on
