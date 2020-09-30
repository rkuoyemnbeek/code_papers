function [U, U_E, Z_E] = lowrank_R_extended( W_U, W_Z, e_R)
% LOWRANK_X_EXTENDED - Gives a low rank representation of the approximate residual
% W_U*W_Z with tolerance e_R 
% 
% Input:
%   (*) W_U: non-orthogonalized basis for the approximate residual
%   (*) W_Z: weights such that W_U*W_Z is the approximate residual
%   (*) e_R: tolerance for the lowrank-approximation
% Output:
%   (*) U: (n x k): orthonormal k-dimensional basis U
%   (*) U_E: part of the basis that is ignored
%   (*) Z_E: weights of the part that is ignored.


%   subfunction of Arnoldi-extended
%   Author: Koen Ruymbeek   
%   Address: Celestijnenlaan 200A, 3001 Leuven
%   email: koen.ruymbeek@cs.kuleuven.be
%   Website: https://www.kuleuven.be/wieiswie/nl/person/00114268
%   Date: 28-Jan-2019; Last revision: 28-Jan-2019
%
%   Copyright (c) 2019, Author
%   All rights reserved.

[Q_U, R_U] = qr( W_U,0);
[Q_Z, R_Z] = qr( W_Z',0);
[U_R, Sigma, Z_R] = svd( R_U*R_Z');
if size(Sigma,2) > 1
    Sigma_diag = diag( Sigma);
    idx = find( cumsum(Sigma_diag(end:-1:1))/sum(Sigma_diag) < e_R, 1, 'last');
    nbr_sigma = length( Sigma_diag);
else
    idx = 0;
    nbr_sigma = 1;
end

if isempty(idx)
    idx = nbr_sigma;
else
    idx = nbr_sigma-idx;
end

U = Q_U*U_R(:,1:idx);
U_E = Q_U*U_R(:,idx+1:end);
Z_E = Sigma(idx+1:end, idx+1:end)*Z_R(:,idx+1:end)'*Q_Z';