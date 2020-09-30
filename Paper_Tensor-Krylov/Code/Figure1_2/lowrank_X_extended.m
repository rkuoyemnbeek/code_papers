function [U,Z, Z_E] = lowrank_X_extended( X_m, V, e_X)
%LOWRANK_X_EXTENDED - Low-rank-approximation of the matrix containing the
%eigenvectors over the parameter. We calculate other vectors such that we
%can do an error analysis.

% Input:
%   (*) X_m: k x n1: [x^V(w1), x^V(w2), ..., x^V(wn1)]: the searched
%       eigenvector of the projected EP for different parameter samples w1,
%       w2,..., n1
%   (*) V: n x k: basis for the k-dimensional subspace 
%   (*) e_X: tolerance for the lowrank-approximation
% Output:
%   (*) U: (k x d1): orthonormal d1-dimensional basis U
%   (*) Z: (d1 x n1): weights such that X_m = U*Z + e_1
%   (*) Z_E: (d2 x n1): associated weights of U_E, with U_E the basis of
%        part that is deleted from the basis.

%   subfunction of Arnoldi-extended
%

%   Author: Koen Ruymbeek   
%   Address: Celestijnenlaan 200A, 3001 Leuven
%   email: koen.ruymbeek@cs.kuleuven.be
%   Website: https://www.kuleuven.be/wieiswie/nl/person/00114268
%   Date: 01-Feb-2019; Last revision: 01-Feb-2019
%
%   Copyright (c) 2019, Author
%   All rights reserved..

n = min( size(X_m));
[Q_X, R_X] = qr(X_m',0);
[U_R, Sigma, Z_R] = svd( R_X');
if size(X_m,2) > 1
    Sigma_diag = diag( Sigma);
    idx = find(cumsum(Sigma_diag(end:-1:1))/sum(Sigma_diag) < e_X, 1, 'last');
else
    idx = 0;
end

if isempty(idx)
    idx = n;
else
    idx = n-idx;
end
U = V*U_R(:,1:idx);
Z = Sigma(1:idx,1:idx)*Z_R(:,1:idx)'*Q_X';

Z_E = Sigma(idx+1:end,idx+1:end)*Z_R(:,idx+1:end)'*Q_X';