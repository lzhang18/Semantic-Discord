function [val] = sampleThreshold(data, context_len, perc, sample_size)
n = length(data)-context_len+1;
idx = randperm(n);
idx1 = idx(1:floor(n/2));
idx2 = idx(floor(n/2)+1:2*floor(n/2));

count = 0;

tmp = abs(idx1-idx2) > context_len;
n_pair = min(sum(tmp), sample_size);


X = zeros(n_pair,context_len);
Y = zeros(n_pair, context_len);
for i = 1: n_pair
    if tmp(i) == 1
        count = count+1;
        X(count,:) = data(idx1(i):idx1(i)+context_len-1);
        Y(count,:) = data(idx2(i):idx2(i)+context_len-1);
    end    
end

X = X(1:count, :);
Y = Y(1:count, :);
diff = zscore(X, 1, 2)- zscore(Y,1,2);
dist = vecnorm(diff,2,2);

dist = sort(dist);
val = prctile(dist,perc*100);

end