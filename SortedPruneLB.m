function [D_BSF,i_BSF,j_BSF,l_BSF, skip_idx, ratio]=SortedPruneLB(data,k,context_len,target_len,context_mu,...
    context_sigma,target_mu,target_sigma, skip_idx, min_sigma, delta)


%n = length(data);

target = data(k:k+target_len-1);


%compute context stats -- O(n)
qt = computeQT(data, target);

% compute lower bound
LB = computeLB2(k,context_len, target_len, target_mu, target_sigma,...
    context_sigma, qt);

LB = real(LB);


% update skip_idx for k
context = data(k:k+context_len-1);
dp_context_vec = MASS_V2(data, context);

% update last and remove first row of V
skip_idx = circshift(skip_idx,[-1 0]);
skip_idx(context_len-target_len+1,:) = (dp_context_vec < delta);

V_colmax = max(skip_idx);   % column max

% set LB of unsimilar contexts to inf for each (k, l)  -- O(n)
V_movmax_l = movmax(V_colmax, [context_len - target_len+1, 0]);
LB(V_movmax_l==0) = Inf;

% compute distance in increasing order of LB
[val,idx]=sort(LB);

D_BSF = Inf;
i_BSF = 0;
j_BSF = 0;
l_BSF = 0;
ratio=1;


if(target_sigma(k)<min_sigma)
    return
end

for l = 1: length(val)
    % skip overlapping target
    if(abs(idx(l)-k)<context_len)
        continue;
    end
    % skip the first context
    if(idx(l)<context_len)
        continue;
    end
    
    % lower bound dist of current iteration
    D_LB=val(l);
    if D_LB > D_BSF
        ratio=l;
        return;
    else
        
        [D, min_i, min_j] = ReconDist_vec(k, idx(l), context_mu, context_sigma, target_mu, target_sigma, qt, context_len, target_len, skip_idx);
        
        if D < D_BSF
            D_BSF = D;
            i_BSF = min_i;
            j_BSF = min_j;
            l_BSF = idx(l);
        end
    end
end
ratio=l;
