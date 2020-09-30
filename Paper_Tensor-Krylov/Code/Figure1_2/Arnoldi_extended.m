function [V, Lambdav, Xv_m, residu_m, E1_m, E2_m, U_v, U12_m, restart_v, nbr_inner_products] = Arnoldi_extended(A_c, norm_A0, f, w_v, U, epsR, epsX, max_rank, max_iter, var_tol_b, eps)
% ARNOLDI_EXTENDED - Extended implementation of the parameteric residual Arnoldi Algoritme 
% for the case A(w) = A0 + f(w)*A1 (1)
% Extra variables containing the error. are calculated and are given back as output
% Figures 1,

% Input:
%   (*) A_c: {A0, A1}: the matrices from (1)
%   (*) f: function from (1)
%   (*) w_v: parameter samples where we want the rightmost eigenvalue from
%   (*) U: n x 1: initial basis for the eigenvectors
%   (*) epsR: the tolerance for the low-rank approximation of the residual
%           (see paper)
%   (*) epsX: the tolerance for the lowrank-approximation of the matrix
%           containing the eigenvectors
%   (*) max_rank: maximal dimension of the subspace.
%   (*) max_iter: maximal nbr of iterations
%   (*) var_tol_b: boolean: 0: fixed eps_X
%                           1: variable eps_X
%   (*) eps: tolerance that we want to achieve

% Output:
%   (*) V: n x k: orthonormal basis for the searched eigenvectors
%   (*) Lambdav: for each sample, the searched eigenvalue
%   (*) Xv_m: for each sample, the eigenvector associated with the searched
%           eigenvalue.
%   (*) residu_m: calculated residual (no low-rank approxmimation)
%   (*) E1_m: for each point and for each iteration, the made error e1
%   (*) E2_m: for eahc point and for each iteration, the made error e2
%   (*) U_v: The size of the subspace over all iterations.
%   (*) U12_m: fore each iteration: (nbr vectors we keep, nbr vectors we
%               throw away)
%   (*) restart_v: vector with nbr of iteration where we neede to restart
%   (*) nbr_inner_products: calculation of the number of inner products
%       that needed to be performed

%   Subfunctions: lowrank_R_extended, lowrank_X_extended
%   MAT-files required: none
%

%   Author: Koen Ruymbeek   
%   Address: Celestijnenlaan 200A, 3001 Leuven
%   email: koen.ruymbeek@cs.kuleuven.be
%   Website: https://www.kuleuven.be/wieiswie/nl/person/00114268
%   Date: 01-Feb-2019; Last revision: 19-maart-2019
%
%   Copyright (c) 2019, Author
%   All rights reserved..


%% Initialize variables
A0 = A_c{1};
A1 = A_c{2};
n = size(A0,1);
n1 = length(w_v);
Lambdav = zeros(n1,1);
residu_m = zeros(n1, max_iter);
count = 0;
U_v = zeros(max_iter,1);
E1_m = zeros(n1, max_iter);
E2_m = zeros(n1, max_iter);
U12_m = zeros(max_iter,2);
hatres_normF = 1;
V = zeros(n, max_rank); 
k = 0;
nbr_iter = 0;
restart_v = zeros(max_iter,1);
count_restart = 1;
restart_v(1) = 0;
f1_wv = f(w_v);
nbr_inner_products = 0;
while any( hatres_normF/sqrt(n1)*norm_A0 > eps) && count < max_iter
    count = count + 1;

    V(:,k+1:k+size(U,2)) = U;

    %% Make reduced eigenvalue problem
    
    if count-1 == restart_v(count_restart)
        pA0_res = A0*V(:,1:size(U,2));
        pA1_res = A1*V(:,1:size(U,2));
        pA0 = V(:,1:size(U,2))'*pA0_res;
        pA1 = V(:,1:size(U,2))'*pA1_res;
        nbr_inner_products = nbr_inner_products + 2*size(U,2)^2;
    else
        pA0_res = [pA0_res, A0*V(:,k+1:k+size(U,2))];
        pA0_help = V(:,1:k+size(U,2))'*pA0_res(:,k+1:k+size(U,2));
        pA0_help1 = [pA0; V(:,k+1:k+size(U,2))'*pA0_res(:,1:k)];
        pA0 = [pA0_help1, pA0_help];

        pA1_res = [pA1_res, A1*V(:,k+1:k+size(U,2))];
        pA1_help = V(:,1:k+size(U,2))'*pA1_res(:,k+1:k+size(U,2));
        pA1_help1 = [pA1; V(:,k+1:k+size(U,2))'*pA1_res(:,1:k)];
        pA1 = [pA1_help1, pA1_help];
        
        nbr_inner_products = nbr_inner_products + 2*( 2*(k+size(U,2))*size(U,2));
    end

    %% Solve reduced eigenvalue problem
    k = k + size(U,2);
    U_v(count) = k;   
    Xv_m = zeros( k, n1);
    Xv_real = zeros(k, n1);
    Xv_comp = Xv_real;
    parfor p = 1:n1
        [xv, lambdav] = eig( pA0+f1_wv(p)*pA1, 'vector');           
        [~, idx] = max(real(lambdav));
        Lambdav(p) = lambdav(idx);
        Xv_m(:,p) = xv(:,idx);
        if norm( imag( xv(:,idx))) > 10^(-13) 
            Xv_real(:,p) = real( xv(:,idx));
            Xv_comp(:,p) = imag( xv(:,idx));                    
        else
            Xv_real(:,p) = Xv_m(:,p);
        end
    end
    complex_v = find( any( Xv_comp ~= 0));
    
    %% Divide in real and complex part
    Xv_all = [Xv_real Xv_comp(:,complex_v)];
    residu_v = vecnorm( pA0_res*Xv_m + pA1_res*Xv_m*diag(f1_wv) - V(:,1:k)*Xv_m*diag(Lambdav)); 
    residu_m(:,count) = residu_v;

    %% low-rank approximation of the eigenvectors of the projected problem
    if k > max_rank && nbr_iter > 2 
        [U,Z, Z_E] = lowrank_X_extended( [Xv_all(:,1:n1), Xv_all(:,n1+(1:length(complex_v)))], V(:,1:k), 10^(-12));
    else
        if var_tol_b == 1
            [U,Z, Z_E] = lowrank_X_extended( [Xv_all(:,1:n1), Xv_all(:,n1+(1:length(complex_v)))], V(:,1:k),  epsX*hatres_normF/(sqrt(n1)*norm_A0));
        else
            [U,Z, Z_E] = lowrank_X_extended( [Xv_all(:,1:n1), Xv_all(:,n1+(1:length(complex_v)))], V(:,1:k), epsX);
        end
    end

    if isempty( complex_v) == 0
        Z(:,complex_v) = Z(:,complex_v) + Z(:,n1+(1:length(complex_v)))*1i; 
        Z_E(:,complex_v) = Z_E(:,complex_v) + Z_E(:,n1+(1:length(complex_v)))*1i;
    end
    Z = Z(:,1:n1);
    %% Determine Error E1    
    Z_E = Z_E(:,1:n1);
    
    E1_m(:,count) = vecnorm( Z_E); 
    U12_m(count,:) = [size(U,2), size(U,2) + size(Z_E,1)];
    
    if k > max_rank && nbr_iter > 2
        k = 0;
        nbr_iter = 0;
        count_restart = count_restart + 1;
        restart_v(count_restart) = count;        
    else
        %% Determine A(X)
        W_U = [A0*U, A1*U];
        W_Z = [Z; Z*diag(f(w_v))];
        %% Orthogonalize X_U
        for i = 1:3
            h = V(:,1:k)'*W_U;
            W_U = W_U - V(:,1:k)*h;
        end   
        hatres_normF = norm(W_U*W_Z, 'fro'); % Frobenius norm of the residu
        
        %% Low-Rank Approximation of X_U X_Z
        norm_W_Z = vecnorm( W_Z);
        W_Z = W_Z*diag(1./norm_W_Z);

        complex_v = find(vecnorm( imag( W_Z)) > 10^(-15));
        W_Z = [real(W_Z), imag(W_Z(:,complex_v))];
        
        [U, U_E, Z_E] = lowrank_R_extended( W_U, W_Z, epsR); 

        if isempty( complex_v) == 0
            Z_E(:,complex_v) = Z_E(:,complex_v) + Z_E(:,n1+(1:length(complex_v)))*1i;
        end
        Z_E = Z_E(:,1:n1);

        E2 = U_E*Z_E;                       
        E2_m(:,count) = vecnorm(E2)';

        for i = 1:3
            h = V(:,1:k)'*U;
            U = U - V(:,1:k)*h;
            U = U*diag(1./vecnorm(U));
        end
        nbr_iter = nbr_iter + 1;
    end

end
E1_m = E1_m(:,1:count);
E2_m = E2_m(:,1:count);
U_v = U_v(1:count);
U12_m = U12_m(1:count,:);
restart_v = restart_v(1:count_restart-1);
residu_m = residu_m(:,1:count);
