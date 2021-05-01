% compute reconstruct lower bound given k
function [LB] = computeLB2(k,context_len, target_len, target_mu, target_sigma,...
    context_sigma, qt)
%compute lower bound
i = k + target_len - context_len;
max_sigma = movmax(context_sigma,[0 context_len - target_len+1],1);

tmp=1:length(max_sigma);
idx=tmp+target_len-context_len;
idx(idx<1)=1;
max_sigma=max_sigma(idx);

max_sigma_i = max_sigma(k);

q = qt./target_len - target_mu(k).*target_mu;
q = q./(target_sigma(k).*target_sigma);

q(q<0)=0;

LB = target_len.*(1-q.^2);
target_sigma = target_sigma(1:length(max_sigma));
LB = LB(1:length(max_sigma));

LB_j = target_sigma./max_sigma.*sqrt(LB);
LB_i = target_sigma(k)/max_sigma_i.*sqrt(LB);
LB_j(LB_j<LB_i)=LB_i(LB_j<LB_i);
LB=LB_j;
end