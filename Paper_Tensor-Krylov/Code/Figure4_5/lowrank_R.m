function [U, Z] = lowrank_R( X_U, X_Z, erel)
% LOWRANK_R - Gives a low rank representation of the matrix X_U*X_V 
% 
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none


%   Author: Koen Ruymbeek   
%   Address: Celestijnenlaan 200A, 3001 Leuven
%   email: koen.ruymbeek@cs.kuleuven.be
%   Website: https://www.kuleuven.be/wieiswie/nl/person/00114268
%   Date: 28-Jan-2019; Last revision: 28-Jan-2019
%
%   Copyright (c) 2019, Author
%   All rights reserved.

[Q_U, R_U] = qr( X_U,0);
[Q_Z, R_Z] = qr( X_Z',0);
[U_R, Sigma, Z_R] = svd( R_U*R_Z');
if size(Sigma,2) > 1
    Sigma_diag = diag( Sigma);
    idx = find( cumsum(Sigma_diag(end:-1:1))/sum(Sigma_diag) < erel, 1, 'last');
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
Z = Sigma(1:idx, 1:idx)*Z_R(:,1:idx)'*Q_Z';

