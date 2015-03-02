
%dat=[x0911(:,2:2:48); x9508_star(:,2:2:48)];
load('dat.mat');
value=dat;
t=isnan(dat);
value(t)=0;
count=dat.*0 + 1;
count(t)=0;

sumofday=sum(value,2);
sumofdayval=sum(count,2);
sod=zeros(365,1);
sodval=zeros(365,1);
for i = 1:27
    sod = sod + sumofday((i-1)*365+1:(i-1)*365+365);
    sodval = sodval + sumofdayval((i-1)*365+1:(i-1)*365+365);
end

aveday = sod./sodval;

d=[31 28 31 30 31 30 31 31 30 31 30 31];
beginday=[1 32 60    91   121   152   182   213   244   274   305   335];
endday=[31 59 90 120 151 181 212 243 273 304 334 365];

som = zeros(1,12);
somval = zeros(1,12);
avemon = zeros(1,12);
for i = 1:12
    som(i) = sum(sod(beginday(i):endday(i)));
    somval(i) = sum(sodval(beginday(i):endday(i)));
    avemon(i) = som(i) / somval(i);
end

% -------------------------------------- Mean Square Fitting ----------------------------------
% compose coefficient matrix
nx=12;
step=30;
A=zeros(365, nx);
for i = 1:nx-1
    for j = 1:step
        prev = (i-1)*step + 1;
        next = i*step + 1;
        x = (i-1)*step + j;
        y1 = (x-prev)/step;
        A(x,i) = y1; 
        A(x,i+1) = 1 - y1;
    end
end
prev = (nx-1)*step + 1;
for x =  prev : 365
    y1 = (x-prev)/(365 - (nx-1)*step);
    A(x,nx) = y1;
    A(x,1) = 1 - y1;
end

Q=A'*A;
[qa qb]=eig(Q);
qt=qa*qb*qa';
fitting=qa*diag(diag(qb).^(-1))*qa'*A'*aveday;
xfit=[1:step:(step*(nx-1)+1) 366];

figure(1);
close Figure 1;
figure(1);
hold on;
plot([1:365], aveday, 'k');
% plot([1:365], aveday2);
% plot([1:365], aveday3);
plot((beginday+endday)/2, avemon, 'r-o','LineWidth',2);
plot(xfit, [fitting; fitting(1)], 'b-*','LineWidth',2);
legend('Solar Radiation', 'Monthly Average Solar Radiation', 'Mean Square Fitting');
xlabel('Day of Year')
ylabel('Average Radiation w/m^2')

% ------------------------------ Average on a certain hour of day --------------

sumofhour=sum(value,1);
sumofhourval=sum(count,1);

avehour = sumofhour./sumofhourval;
figure(2);
close Figure 2;
figure(2);
hold on;
plot([1:24], avehour);
beginy = 25;
stepy = 1;
endy = 27;
% for i=1:length(value(:,1))
%     delta=(i - length(value(:,1))/2) /3000;
%     plot(delta+[1:24],value(i,:),'.', 'MarkerSize',2);
% end
begind0=(beginy-1)*365+1;
endd0=(endy)*365;
for y = beginy:stepy:endy;
    begind = (y-1)*365+1;
    endd = (y+stepy-1)*365;
    for i=begind:endd
        delta = ((i - begind0) / (endd0 - begind0) - 1/2) / 1.5;
        plot(delta+[1:24],value(i,:),'.', 'MarkerSize',2);
    end
end
xlabel('Time of Day')
ylabel('Average Radiation w/m^2')


