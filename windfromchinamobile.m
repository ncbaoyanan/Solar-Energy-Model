clear all
close all
clc

load wind.mat
load wind1.mat

data  = [wind1;wind];
figure;
plot(data,'.-');
xlabel('hour');
ylabel('Wind Strength (Voltage)');
grid on;

[h,center] = hist(data,5);
h = h/sum(h);
figure;
bar(center,h);
xlabel('Wind Strength (Voltage)');
ylabel('PDF');
grid on;

figure;
data1 = data - mean(data);
[c_ww,lags] = xcorr(data1,'coeff');
stem(lags,c_ww)
xlabel('hour');
ylabel('Correlation');
grid on;