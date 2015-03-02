clear all
close all
clc
O = 10;
Q = 10;

load hourno0.mat
data = hourno0;
scale = (max(data)-min(data))/O;
len = length(data);
mindata = min(data);
for i = 1:len
    data(i) = max(ceil((data(i)-mindata)/scale),1);
end
data = data';
% data = reshape(data,365,27);


loglikopt = -inf;
num =  100;
for i = 1:num
% initial guess of parameters
prior1 = normalise(rand(Q,1));
transmat1 = mk_stochastic(rand(Q,Q));
obsmat1 = mk_stochastic(rand(Q,O));

% improve guess of parameters using EM
[LL, prior2, transmat2, obsmat2] = dhmm_em(data, prior1, transmat1, obsmat1,'thresh', 1e-8, 'max_iter', 20);
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