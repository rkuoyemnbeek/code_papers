function [V, X_t, lambda, residu_v] = PSIRA_dD( A_c, F_t, U, epsR, epsX, max_rank, max_iter, eps)
%PSIRA_dD - Parametric Shift-and-Invert Arnoldi:
% Parametrized eigenvalue problem of the form
% A(w) = \sum_{i=0}^m f_i(w) A_i with:
%   A_0: full rank
%   A_i, i > 0: rank-1 matrices
%   f_i, i= 0, \hdots, m: \Omega --> R functions
% Aim: find a basis V such that for all (w1, ..., wd) \in \Omega, the
%   eigenvalue with real part closest to zero is found with certain
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

%% Initialise all variables
n = size(A_c{1},1);
d = length( size( F_t))-1;
m = length(A_c)-1;
count = 0;
norm_residu = 1;
nbr_iter = 0;
F_tt = tt_tensor( F_t);

n_v = zeros(d,1);
u_m = zeros(n,m);
q_m = zeros(n,m);
abs_tol_v = zeros( max_iter,1);
residu_v = ones( max_iter,1);
V = zeros(n, max_rank);

for i = 1:d
    n_v(i) = size( F_tt, i+1);
end
p = prod(n_v);
n_c = num2cell(n_v);
lambda_schat_t = zeros( max_iter, n_c{:});
if d == 1
   lambda_schat_t = zeros( max_iter, n_c{1});
end

%% decompose the rank-1 matrices A_i, i > 0
for i = 1:m
    [u, q]= qr( A_c{i+1});
    j = find( abs( sum(q,2)) > 0);
    u_m(:,i) = u(:,j);
    q_m(:,i) = q(j,:)';
end
A0_u_m = orth( A_c{1}\u_m);
m = size( A0_u_m, 2);

%% orthogonalise start subspace U w.r.t V
V(:,1:m) = A0_u_m;
for i = 1:3
    h = V(:,1:m)'*U;
    U = U - V(:,1:m)*h;
    U = U*diag(1./vecnorm(U));
end
U = orth( U);
k = m;

norm_A0 = normest(A_c{1});
restart_b = 1;
while norm_residu > eps*sqrt(p)*norm_A0 && count < max_iter
    count = count + 1;
    size_U = size(U,2);
    V(:,k+1:k+size_U) = U;
    k = k+size_U;


    %% Solve projected eigenvalue problem
    pA_help_c = cellfun( @(x) x*V(:,1:k), A_c, 'UniformOutput',false);
    pA_c = cellfun( @(x) V(:,1:k)'*x, pA_help_c, 'UniformOutput',false);
    
    X_t = zeros(k,n_c{:}); % the first eigenvector is the vector a that we use (see section 5.1) 
    f1 = F_t(:,1);
    Ahelp = cellfun( @(x,y) x*y, num2cell(f1)', pA_c, 'Uniform', false);
    A = sum( cat(3, Ahelp{:}),3);
    [X, Lambda] = eig( A, 'vector');
    X = X(:,imag( Lambda) >= 0);
    Lambda = Lambda( imag(Lambda) >= 0); 
    % To detect eigenvalues of the reduced problem which accidently came as first in our criterium
    %  we save the first 5 eigenvectors for every sample point, 
    % More info in the paper at section 6.1
    if restart_b == 1
        X_old_c = cell(4,1);
        X_old_c{1} = X_t;
        X_old_c{2} = X_t;
        X_old_c{3} = X_t;
        X_old_c{4} = X_t;
        [~, idx] = sort(abs(real(Lambda)), 'ascend');
        lambda_schat_t(count, 1) = Lambda(idx(1));
        x1 = X(:,idx(1));
        X_t(:,1) = X(:,idx(1));
        for i = 1:min(4:length(Lambda)-1)
            X_old_c{i}(:,1) = X(:,i+1);
        end
    else
        X_c = cell(4,1);
        X_c{1} = X_t;
        X_c{2} = X_t;
        X_c{3} = X_t;
        X_c{4} = X_t;
        [~, idx_sort] = sort(abs(real(Lambda)), 'ascend');
        i = 1;
        while i < length( Lambda)+1
            if  max( abs(   X(1:k-size_U,idx_sort(i))'*([X_t_old(:,1), X_old_c{1}(:,1), X_old_c{2}(:,1), X_old_c{3}(:,1), X_old_c{4}(:,1)] ))) < 0.99
                i = i+1;
            else
                break
            end
        end
        if i == length(Lambda)+1
            i = 1;
        end
        
        for j = 1:min( i-1,4)
            X_c{j}(:,1) = X(:,idx_sort(j));
        end
        for j = (i+1):min(5, length(Lambda))
            X_c{j-1}(:,1) = X(:,idx_sort(j));
        end            
        lambda_schat_t(count, 1) = Lambda(idx_sort(i));
        x1 = X(:,idx_sort(i));
        X_t(:,1) = X(:,idx_sort(i));        
    end
    
    for p1 = 2:prod(n_v)
        f1 = F_t(:,p1);
        Ahelp = cellfun( @(x,y) x*y, num2cell(f1)', pA_c, 'Uniform', false);
        A = sum( cat(3, Ahelp{:}),3);
        [X, Lambda] = eig( A, 'vector'); 
        X = X(:,imag( Lambda) >= 0);
        Lambda = Lambda( imag(Lambda) >= 0);
        if restart_b == 1

            [~, idx] = sort(abs(real(Lambda)), 'ascend');
            lambda_schat_t(count, p1) = Lambda(idx(1));
            X_t(:,p1) = X(:,idx(1));
            for i = 1:min(4:length(Lambda)-1)
                X_old_c{i}(:,p1) = X(:,i+1);
            end
        else
            [~, idx_sort] = sort(abs(real(Lambda)), 'ascend');
            i = 1;
            while i < length( Lambda)+1
                if  max( abs(   X(1:k-size_U,idx_sort(i))'*([X_t_old(:,p1), X_old_c{1}(:,p1), X_old_c{2}(:,p1), X_old_c{3}(:,p1), X_old_c{4}(:,p1)] ))) < 0.99
                    i = i+1;
                else
                    break
                end
            end
            if i == length(Lambda)+1
                i = 1;
            end

            for j = 1:min( i-1,4)
                X_c{j}(:,p1) = X(:,idx_sort(j));
            end
            for j = (i+1):min(5, length(Lambda))
                X_c{j-1}(:,p1) = X(:,idx_sort(j));
            end            
            lambda_schat_t(count, p1) = Lambda(idx_sort(i));
            X_t(:,p1) = ( X(:,idx_sort(i))'*x1)/abs(X(:,idx_sort(i))'*x1)*X(:,idx_sort(i));       
        end
    end  
    
    X_t_old = X_t;
    if restart_b == 1
        restart_b =0;
    else
        X_old_c = X_c;
    end
   %% 

    if k > max_rank && nbr_iter > 2
        %% low-rank approximation of the eigenvectors over the parameter
        abs_tol_v(count) = epsX*norm_residu/sqrt(p)/norm_A0;
        U = tt_tensor_split( X_t, abs_tol_v(count)); 
        residu_v(count) = residu_v(count-1);

        U = V(:,1:k)*U;
        for i = 1:3
            h = V(:,1:m)'*U;
            U = U - V(:,1:m)*h;
            U = U*diag(1./vecnorm(U));
        end
        U = orth(U);
        k = m;

        nbr_iter = 0;
        restart_b = 1;
    else
        %% low-rank approximation of the eigenvectors over the parameter
        if count == 1
            abs_tol_v(count) = epsX;
        else
            abs_tol_v(count) = min( [epsX, epsX*norm_residu/sqrt(p)]);
        end
        [U, Z_tt] = tt_tensor_split( X_t, abs_tol_v(count));
        U = V(:,1:k)*U;
        
        %% Determine the residual 
        W_U = [A_c{1}*U, u_m];
        W_Z_tt = W_Z_maker_rank1( Z_tt, F_tt, U, q_m);
        W_Z_tt = round( W_Z_tt, 10^(-12));

        for i = 1:3 
           h = V(:,1:k)'*W_U;
           W_U = W_U-V(:,1:k)*h;
        end 
        residu_tt = ttm( W_Z_tt, 1, W_U');

        residu_v(count) = norm( residu_tt);      
        norm_residu = residu_v(count);
        %% low-rank approximation of the residual
        U =round_rightleft_1mode(W_U, W_Z_tt, epsR);
        
        %% solving system
        U_help = A_c{1}\U;      
        for i = 1:3
            h = V(:,1:k)'*U_help;
            U_help = U_help - V(:,1:k)*h;
            U_help = U_help*diag(1./vecnorm(U_help));
        end
        U = orth(U_help);  
        nbr_iter = nbr_iter + 1;
    end  
end
if count < max_iter
    display( 'it is converged!')
else
    display( 'not converged')
end

residu_v = residu_v(1:count);

U= tt_tensor_split( X_t, 10^(-12));
V = V(:,1:k)*U;
lambda = lambda_schat_t(count,:);