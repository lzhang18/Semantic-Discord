function [D_BSF,i_BSF,j_BSF,l_BSF, ratio, delta] = SemanticDiscord(data, context_len, target_len,min_sigma)
%main function for semantic discord. The function requires following inputs
%
% data: 1-D time series, size Nx1
%
% context_len: context length (see paper for detail)
%
% target_len: target length (see paper for detail)
%
% min_sigma: a variance threshold. Any subsequece has
%            variance smaller than min_sigma will be
%            ignore. This is very useful in practical.

%remove NaN
data(isnan(data))=0;

%Compute distance threshold based on random sampled subsequence
%the threshold will be used to filter out dissimilar context

sample_size = 2000; % sample size
threshold = 0.4;

% compute distance threshold
delta =  sampleThreshold(data, context_len, threshold, sample_size);


% compute context and target mean and standard dev
context_mu = movmean(data,[context_len-1 0]);
context_mu = context_mu(context_len:end);
context_sigma = movstd(data,[context_len-1 0],1);
context_sigma = context_sigma(context_len:end);
target_mu = movmean(data,[target_len-1 0]);
target_mu = target_mu(target_len:end);
target_sigma = movstd(data,[target_len-1 0],1);
target_sigma = target_sigma(target_len:end);

% compute skip_idx, where any subseuences has distance greater than delta
% is 0

skip_idx = zeros(context_len-target_len+1, length(data)-context_len+1);
dp_context = zeros(context_len-target_len+1, length(data)-context_len+1);

k = context_len;
for c = 1: context_len-target_len+1
    context_start = k+target_len-context_len+c-1;
    context = data(context_start:context_start+context_len-1);
    dp_context_vec = MASS_V2(data, context);
    dp_context(c,:) = dp_context_vec;
    skip_idx(c, :) = (dp_context_vec < delta);
end

%start main function
tic

for inc=1:length(data)-2*context_len
    k=inc+context_len;
    [D_BSF(k),i_BSF(k),j_BSF(k),l_BSF(k),skip_idx, ratio(k)]=SortedPruneLB(data,k,context_len,target_len,context_mu,context_sigma,target_mu,target_sigma, skip_idx, min_sigma, delta);
    disp(k) 
end
toc