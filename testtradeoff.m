clear all
close all
clc
f := stats::erlangCDF(2, 1):
f(-infinity), f(-3), f(0.5), f(2/3), f(PI), f(infinity)
r = 5;
lambda = 1;
% f = @(t) 1 - exp(-r*lambda*t)*sum((r*lambda*t).^(0:r-1)./factorial(1:r-1));

Nt = 5;
Tset = exp(-1:0.1:1);
Pi = f(Tset);
Ei = f(Tset);
Len = length(Tset);
for i1 = 1:Len
    for i2 = 1:Len;
        for i3 = 1:Len
            for i4 = 1:Len
                for i5 = 1:Len
                end
            end
        end
    end
end