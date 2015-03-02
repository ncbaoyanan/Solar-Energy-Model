clear all
close all
clc
O = 5;
Q = 5;

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

data = seg;


loglikopt = -inf;
num =  1;
for i = 1:num
% initial guess of parameters
prior1 = normalise(rand(Q,1));
transmat1 = mk_stochastic(rand(Q,Q));
obsmat1 = mk_stochastic(rand(Q,O));

% improve guess of parameters using EM
[LL, prior2, transmat2, obsmat2] = dhmm_em(data, prior1, transmat1, obsmat1,'thresh', 1e-8, 'max_iter', 50);
% LL

% use model to compute log likelihood
loglik = dhmm_logprob(data, prior2, transmat2, obsmat2);
% log lik is slightly different than LL(end), since it is computed after the final M step
    if loglik > loglikopt
        prioropt = prior2;
        transmatopt = transmat2;
        obsmatopt = obsmat2;
        loglikopt = loglik;
    end
    i/num
end

loglikopt
prioropt
transmatopt
obsmatopt

figure;
mesh(obsmatopt)
title('Observation matrix')
xlabel('State now')
ylabel('Observation now')

figure;
mesh(transmatopt)
title('Transition matrix')
xlabel('State now')
ylabel('State next')