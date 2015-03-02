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
    end
    i = i + 1;
end

seg = seg(1:Nseq);

data = seg;

num =  100;
cores = 10;

Mloglikopt = cell(cores,1);
Mprioropt = cell(cores,1);
Mtransmatopt = cell(cores,1);
Mobsmatopt = cell(cores,1);
for k = 1:cores
    kk  = 0;
    Mloglikopt(k) = {-inf};
    for kk = 1:num
    % initial guess of parameters
    prior1 = normalise(rand(Q,1));
    transmat1 = mk_stochastic(rand(Q,Q));
    obsmat1 = mk_stochastic(rand(Q,O));

    % improve guess of parameters using EM
    [LL, prior2, transmat2, obsmat2] = dhmm_em(data, prior1, transmat1, obsmat1,'thresh', 1e-8, 'max_iter', 20,'adj_prior',0,'verbose', 0);
    % LL

    % use model to compute log likelihood
    loglik = dhmm_logprob(data, prior2, transmat2, obsmat2);
    % log lik is slightly different than LL(end), since it is computed after the final M step
        if loglik > Mloglikopt(k)
            Mprioropt(k) = {prior2};
            Mtransmatopt(k) = {transmat2};
            Mobsmatopt(k)  = {obsmat2};
            Mloglikopt(k) = {loglik};
        end
    end    
end


[loglikopt, index] = min(Mloglikopt);
prioropt = Mprioropt(index);
transmatopt = Mtransmatopt(index);
obsmatopt = Mobsmatopt(index);

figure;
surf(obsmatopt)
title('Observation matrix')
xlabel('State now')
ylabel('Observation now')

figure;
surf(transmatopt)
title('Transition matrix')
xlabel('State now')
ylabel('State next')