function [U,Z] = lowrank_X( X_m, erel)
%LOWRANK_X - Lowrank-presentation of the matrix containing the eigenvectors
% for different parameter values. U is a basis for the column space. 

%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%
%   See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

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
    idx = find( cumsum(Sigma_diag(end:-1:1))/sum(Sigma_diag) < erel, 1, 'last');
else
    idx = 0;
end

if isempty(idx)
    idx = n;
else
    idx = n-idx;
end
U = U_R(:,1:idx);
Z = Sigma(1:idx,1:idx)*Z_R(:,1:idx)'*Q_X';
end
