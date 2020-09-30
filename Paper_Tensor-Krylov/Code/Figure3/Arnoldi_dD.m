function [V, X_schat_t, lambda_schat_t, residu_v] = Arnoldi_dD( A_c, F_t, F_tt, w_c, U, tol_W, tol_X, max_rank, eps)
%ARNOLDI_dD - Parametric Residual Arnoldi for multiple parameters:
% Parametrized eigenvalue problem of the form
% A(w) = \sum_{i=0}^m f_i(w) A_i with:
%   A_0: full rank
%   A_i, i > 0: rank-1 matrices
%   f_i, i= 0, \hdots, m: \Omega --> R functions
% Aim: find a basis V such that for all (w1, ..., wd) \in \Omega, the
%   eigenvalue with largest real part is found with certain
%   accuracy. Let \Omega \subset \R^d, and we first discretize each dimension 
%   i in n_i elements {w^i_1, \hdots, w^i_{n_i} }. We then make a
%   d-dimensional grid using kronecker product, such that we have in total
%   n_1 n_2 ... n_d parameter values. For each element we approximate 
%   the eigenvalue with real part closest to zero.
%   

% INPUT:
%   (*) A_c: cell(m+1,1) = {A_0, A_1, \hdots, A_m}
%   (*) F_t: (m+1) x n_1 x \hdots x n_d real tensor with:
%           F_t(i0,i1,\hdots,id) = f_i0( w^1_i1, w^2_i2, \hdots, w^d_id)   
%   (*) U: n x k: orthonormale startbasis
%   (*) epsR: tolerance for the low-rank approx of the residual over the
%            parameter
%   (*) epsX: tolerance for the low-rank approx of the eigenvector over the
%               parameter
%   (*) max_rank: maximal rank of the subspace (not higher than 100)
%   (*) max_iter: maximal nbr of iterations
%   (*) eps = not less than 10^(-10), the higher it is, the less iterations
%   need to be done
% OUTPUT:
%   (*) V: n x p: the resulting basis, can be use 
%   (*) X_t = p x n1 x \hdots x nd: tensor with eigenvectors of the projected problem
%           V^T A(w) V for each w in grid. In each X_t(:,i1,...,id) the searched
%           eigenvector of V^T A(w^1_i1, w^2_i2, \hdots, w^d_id) V is placed. The
%            approximate eigenvector of A(w^1_i1, w^2_i2, \hdots, w^d_id) is V*X_t(:,i1,...,id)
%   (*) lambda: n1 x n2 x \hdots x nd: for each parameter value on the
%       grid, an estimation of the searched eigenvalue.
%   (*) residu_v: vector with for all iterations, the norm of the residual
%       over all parameter values, each element is 
%       \sqrt( sum_{all parameter values} r(w^1_i1, \hdots, w^d_id)^2)
%       with r(w^1_i1, \hdots, w^d_id) the residual for parameter value 
%       (w^1_i1, \hdots, w^d_id)

%   packages required: TT-Toolbox (= Tensor-Train Toolbox by Oseledets)
%   Subfunctions: tt_tensor_split, round_rightleft_1mode,
%       W_Z_maker_rank1, ttm (in TT-toolbox)
%   MAT-files required: none

%   Author: Koen Ruymbeek   
%   Address: Celestijnenlaan 200A, 3001 Leuven
%   email: koen.ruymbeek@cs.kuleuven.be
%   Website: https://www.kuleuven.be/wieiswie/nl/person/00114268
%   Date: 24-Sep-2019; Last revision: 24-Sep-2019
%
%   Copyright (c) 2019, Author
%   All rights reserved.

d = length(w_c);
n_v = zeros(d,1);
for i = 1:d
    n_v(i) = length(w_c{i});
end
n = size(A_c{1},1);
p = prod(n_v);
n_c = num2cell(n_v);
max_iter = 1000;
residu_v = ones( max_iter,1);

teller_herstart = 1;

k = 0;
nbr_iter = 0;

lambda_schat_t = zeros( max_iter, n_c{:});
if d == 1
    lambda_schat_t = zeros( max_iter, n_c{1});
end
teller = 0;
norm_residu = 1;
norm_A0 = normest(A_c{1});
V = zeros(n, max_rank);

while norm_residu > eps*sqrt(p)*norm_A0 && teller < max_iter
    teller = teller + 1

    V(:,k+1:k+size(U,2)) = U;
    [k size(U,2)]
    k = k+size(U,2);

    %% geprojecteerd eigenw probleem oplossen  
    G_res_c = cellfun( @(x) x*V(:,1:k), A_c, 'UniformOutput',false);
    G_c = cellfun( @(x) V(:,1:k)'*x, G_res_c, 'UniformOutput',false);

    X_schat_t = zeros(k,n_c{:});
    ones_v = ones( k,1);
    
    parfor p1 = 1:prod(n_v)
        f1 = F_t(:,p1);
        Ahelp = cellfun( @(x,y) x*y, num2cell(f1)', G_c, 'Uniform', false);
        A = sum( cat(3, Ahelp{:}),3);
        [X, Lambda] = eig( A, 'vector'); 
        [~, idx] = min(real(Lambda));
        lambda_schat_t(teller, p1) = Lambda(idx);
        X_schat_t(:,p1) = sign( X(:,idx)'*ones_v)*X(:,idx);
    end 

    if k > max_rank && nbr_iter > 2
        [U, Z_tt] = tt_tensor_split( X_schat_t, min( [tol_X, tol_X*norm_residu/sqrt(p)]));
        U = V(:,1:k)*U; 
        residu_v(teller) = residu_v(teller-1);
    else
        [U, Z_tt] = tt_tensor_split( X_schat_t, min( [tol_X, tol_X*norm_residu/sqrt(p)])); U = V(:,1:k)*U;     
    end
    if k > max_rank && nbr_iter > 2      
        nbr_iter = 0;
        teller_herstart = teller_herstart + 1;
        k = 0;
    else
        W_U = cell2mat( cellfun( @(y) y*U, A_c, 'UniformOutput',false));
        W_Z_tt = W_Z_maker(Z_tt, F_tt);
        W_Z_tt = round( W_Z_tt, 10^(-12));
        for i = 1:3
            h = V(:,1:k)'*W_U;
            W_U = W_U-V(:,1:k)*h;
        end        
        residu_tt = ttm( W_Z_tt, 1, W_U');
        residu_v(teller) = norm( residu_tt); 
        
        norm_residu = residu_v(teller);
        U=round_rightleft_1mode(W_U, W_Z_tt, tol_W);          
        for i = 1:3
            hU = V(:,1:k)'*U;
            U = U-V(:,1:k)*hU;
            U = U*diag(1./vecnorm(U));
        end
        nbr_iter = nbr_iter + 1;
    end   
end

[U, ~] = tt_tensor_split( X_schat_t, min( [tol_X*norm_residu/sqrt(p),10^(-12)]));
V = V(:,1:k)*U; 
residu_v = residu_v(1:teller);
