% vectorized version of reconstruct distance matrix given k, l

function [d, min_i, min_j] = ReconDist_vec(k, l, context_mu, context_sigma, target_mu, target_sigma, ...
    qt, context_len, target_len, V)

i_vec = k+ target_len -context_len: k;
j_vec = l+target_len-context_len:l;
reconstruct_dist_mtx = ReconDistGivenIJ_vec(i_vec,j_vec,k,l,context_mu,...
            context_sigma,target_mu, target_sigma, qt, target_len);

% set dist = inf for unsimilar context pairs
tmp = (V(:, j_vec) == 0);
reconstruct_dist_mtx(tmp) = inf;

% compute minimum distance -- O(context_len^2)
[M, I] = min(reconstruct_dist_mtx);
[d, col_idx] = min(M);

row_idx = I(col_idx);
min_i = k + target_len - context_len+ row_idx-1;
min_j = l + target_len - context_len+ col_idx-1;


end