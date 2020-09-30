function [V, Lambda] = calcul_eigv(A,B,sm_big, stan_gen, opts)
% Calculates the <sm_big> eigenpairs of the matrix-pair (A,B)
% INPUT:
%   (*) (A, B): matrix-pair, if stan_gen = stan, then only A is used
%   (*) sm_big: 'sm' = calculate the smallest eigenvalues, 
%               'big' = calculate the largest eigenvalues
%   (*) stan_gen:   'stan': standard eigenvalue problem is solved
%                   'gen': generalised eigenvalue problem is solved
%   (*) opts: options who are used as options for 'eigs'
% OUTPUT:
%   (*) V: Krylov-space used for calculating the eigenvalue
%   (*) Lambda: Calculated eigenvalue

if strcmp(sm_big, 'sm')
    if strcmp( stan_gen, 'gen')
        [V, Lambda, ~, ~, flag] = eigs_adapt(A, B, 1, 'sa', opts);
        if flag == 1               
            warning('eigs heeft niet geconvergeerd')
            [V, Lambda] = eig(full(A),full(B));
            [Lambda, idx] = sort( diag(Lambda));
            V = V(:,idx(1:20));
            Lambda = Lambda(1);
        end
    else
        
        [V, Lambda,  ~, ~, flag] = eigs_adapt(A, 1, 'sa', opts);
        if flag == 1
            warning('eigs heeft niet geconvergeerd')
            [V, Lambda] = eig(full(A));
            [Lambda, idx] = sort( diag(Lambda));
            V = V(:,idx(1:20));
            Lambda = Lambda(1);
        end
    end

else
    if strcmp( stan_gen, 'gen')
        [V, Lambda,  ~, ~, flag] = eigs_adapt(A, B,1, 'la', opts);
        if flag == 1
            [V, Lambda] = eig(full(A),full(B));
            [Lambda, idx] = sort( diag(Lambda), 'descend');
            V = V(:,idx(1:20));
            Lambda = Lambda(1);
        end
    else
        A = (A+A')/2;
        [V, Lambda,  ~, ~, flag] = eigs_adapt(A, 1, 'la', opts);
        if flag == 1
            [V, Lambda] = eig(full(A));
            [Lambda, idx] = sort( diag(Lambda), 'descend');
            V = V(:,idx(1:20));
            Lambda = Lambda(1);
        end
    end
end