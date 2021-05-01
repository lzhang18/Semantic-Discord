%compute inner products of targets O(n log n)
function [qt] = computeQT(data,y)
target_len = length(y);
n = length(data);
y = y(end:-1:1);%Reverse the query
y(target_len+1:n) = 0; %aappend zeros

%The main trick of getting dot products in O(n log n) time
X = fft(data);
Y = fft(y);
Z = X.*Y;
qt = ifft(Z);
qt = qt(target_len:end);
end