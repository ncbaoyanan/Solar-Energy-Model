clear all
close all
clc
load('dat.mat');
value=dat;

[lx ly] = size(value);

Eseq = zeros(lx*ly,1);

for i = 1:lx
    for j = 1:ly
        Eseq((i-1)*24+j) = value(i,j);
        if ~( Eseq((i-1)*24+j) >= 0 && Eseq((i-1)*24+j) < 1200)
            Eseq((i-1)*24+j) = 0;
        end
    end
end

Year27 = Eseq;
figure;
plot(Year27,'.')
xlabel('Hour')
ylabel('E(t)')

Days = 5;
Len27 = 365*24*27;
Year27D = zeros(Len27,2*Days+1);
for i = 1:2*Days+1
    ii = i - Days;
    Year27D(:,i) = [Year27(mod(1+24*ii,Len27):Len27); Year27(1:mod(24*ii,Len27))];
end

Year27DD = reshape(Year27D,Len27*(2*Days+1),1);
Year27DD = reshape(Year27DD,365*24,27*(2*Days+1));
MaxD = zeros(365*24,1);
ll = 10;
for i= 1:365*24
    s = sort(Year27DD(i,:),'descend');
    MaxD(i) = s(ll);
end
figure;
plot(MaxD,'.r')

Maxhour = reshape(MaxD,24,365);
Maxhour = Maxhour';
figure;
plot(Maxhour);



