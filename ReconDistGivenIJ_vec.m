function [d] = ReconDistGivenIJ_vec(i_vec,j_vec,k,l,context_mu,context_sigma,target_mu, target_sigma,...
    qt, target_len)
%compute context-aware distance
no_of_context = length(i_vec);

context_mu_i_vec = repmat(context_mu(i_vec)',no_of_context,1)';
context_mu_j_vec = repmat(context_mu(j_vec),1, no_of_context)';
context_sigma_i_vec = repmat(context_sigma(i_vec)',no_of_context,1)';
context_sigma_j_vec = repmat(context_sigma(j_vec),1, no_of_context)';

d1_tmp = target_sigma(k)^2+ (target_mu(k)-context_mu(i_vec)).^2;
d1 = target_len./(context_sigma(i_vec).*context_sigma(i_vec)).*d1_tmp;


d2_tmp = qt(l)/target_len;
d2_tmp = d2_tmp - context_mu_i_vec*target_mu(l)- context_mu_j_vec.*target_mu(k)+ ...
    context_mu_i_vec.*context_mu_j_vec;
d2_tmp1 =  -2*target_len./(context_sigma(i_vec).*context_sigma(j_vec)');
d2_mtx = d2_tmp1.*d2_tmp;

d3_tmp = target_sigma(l)^2+ (target_mu(l)-context_mu(j_vec)).^2;
d3 = target_len./(context_sigma(j_vec).*context_sigma(j_vec)).*d3_tmp;

d1_mtx = repmat(d1, 1,no_of_context);
d3_mtx = repmat(d3', no_of_context,1);

d = d1_mtx+d2_mtx+d3_mtx;
d = sqrt(d);
end