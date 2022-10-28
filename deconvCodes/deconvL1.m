function [x, cost] = deconvL1(y, lam, b, a, Nit)
% https://eeweb.engineering.nyu.edu/iselesni/lecture_notes/sparse_deconv/index.html
% [x, cost] = deconvL1(y, lam, b, a, Nit)
% Sparse deconvolution by L1 norm minimization
% Cost function : 0.5 * sum(abs(y-H(x)).^2) + lam * sum(abs(x));
%
% INPUT
%   y - noisy data
%   b, a - filter coefficients for LTI convolution system H
%   lam - regularization parameter
%   Nit - number of iterations
%
% OUTPUT
%   x - deconvolved sparse signal
%   cost - cost function history

% Ivan Selesnick, selesi@poly.edu, 2012
% Algorithm: majorization-minimization with banded system solver.


y = y(:);                                               % convert column vector
cost = zeros(1, Nit);                                   % cost function history
N = length(y);

Nb = length(b);
Na = length(a);
B = spdiags(b(ones(N,1), :), 0:-1:1-Nb, N, N);          % create sparse matrices
A = spdiags(a(ones(N,1), :), 0:-1:1-Na, N, N);

H = @(x) A\(B*x);                                       % filter
AAT = A*A';                                             % A*A' : sparse matrix
x = y;                                                  % initialization
g = B'*(A'\y);                                          % H'*y

for k = 1:Nit   
    W = (1/lam) * spdiags(abs(x), 0, N, N);                 % W : diag(abs(x)) (sparse)
    F = AAT + B*W*B';                                       % F : banded matrix (sparse)
    x = W * (g - (B'*(F\(B*(W*g)))));                       % update
    cost(k) = 0.5*sum(abs(y-H(x)).^2) + lam*sum(abs(x));    % cost function value
end

cost = cost';

