clear all
close all
clc
O = 10;
Q = 10;

scale = 1/O;

load('hour.mat');
len  = length(hour);
for i = 1:len
    if hour(i) > 1
        hour(i) = 0;
    end
end

seg = cell(len,1);

Nseq = 0;
i = 1;
while i < len
    k = 0;
    while hour(i)
        k = k + 1;
        i = i + 1;
        if i > len
            break;
        end
    end
    if k
        Nseq = Nseq + 1;
        seg(Nseq) = {ceil(hour(i-k:1:i-1)/scale)};
        %seg(Nseq) = {hour(i-k:1:i-1)};
    end
    i = i + 1;
end

seg = seg(1:Nseq);

figure;
hold on;
for i = 1:1000
    plot(cell2mat(seg(i)),'.-')
end
xlabel('Time');
ylabel('State')

data = seg;

Trans = zeros(O,O);
for i  = 1:Nseq
    one = cell2mat(seg(i));
    len = length(one);
    while len > 1
        Trans(one(len-1),one(len)) = Trans(one(len-1),one(len)) + 1;
        len = len - 1;
    end    
end

for i = 1:O
    ss = 0;
    for j = 1:O
        ss = ss + Trans(i,j);
    end
    Trans(i,:) = Trans(i,:)/ss;
end

figure;
mesh(Trans)
title('Transition matrix')
xlabel('State now')
ylabel('State next')

