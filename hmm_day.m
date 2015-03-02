clear all
close all
clc
O = 5;
Q = 5;

load day
scale = (max(day)-min(day))/O;
len = length(day);
data = zeros(len,1);
for i = 1:len
    data(i) = max(ceil((day(i)-min(day))/scale),1);
end
data = data';
% data = reshape(data,365,27);


loglikopt = -inf;

for i = 1:10
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
        logliopt = loglik;
    end
end

logliopt
prioropt
transmatopt
obsmatopt