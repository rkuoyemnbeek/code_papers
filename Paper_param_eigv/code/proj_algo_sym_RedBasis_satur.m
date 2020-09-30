function [basis, omega_added, lambda1_m, timings_deriv, timings_eigenv] = proj_algo_sym_RedBasis_satur(matrices, functions, max_vect, nbr_vect, startpoints_c, nbr_valid, tol, sm_big, stan_gen, method, tol_gmres, opts_eigs)
% Aim: Calculating a subspace V such that | \lambda_1^V - \lambda_1| < tol
% INPUT:
%   (*) (matrices, f_c): ({A_c, B_c}, f_c{i} matrices such that A(w) = sum f_c{i}(w) A_c{i}
%           (same for B)
%   (*) max_vect: maximal dimension of the subspace
%   (*) nbr_vect: nbr of (approx) eigenvectors that is added in each
%           iteration
%   (*) startpoints_c: {leftbound, rightbound, nbr_discr_points}
%   (*) nbr_valid: nbr of validation points
%   (*) tol: tolerance level
%   (*) sm_big: 'sm' = calculate the smallest eigenvalues, 
%               'big' = calculate the largest eigenvalues 
%   (*) stan_gen:   'stan': standard eigenvalue problem is solved
%                   'gen': generalised eigenvalue problem is solved
%   (*) method: method for calculating the derivatives: 
%               'nonderiv': no derivatives added
%               'system': derivatives calculated by solving a system
%   (*) tol_gmres: tolerance for the GMRES-method (if we use it)
%   (*) opts_eigs: opts used in eigs_adapt
% OUTPUT:
%   (*) basis: orthonormal basis for the searched subspace
%   (*) omega_added: the parameter values from which the subspace is made
%               (start parameter values included)
%   (*) lambda1_m: eigenvalues for all parametervalues
%   (*) timings_deriv: timings when calculating the derivatives of the
%                       eigenvectors
%   (*) timings_eigenv: timings when calculating the eigenvectors
%   

% zelde als proj_algo_sym_RedBasis maar nu met de saturatie-idx.
% Dit is gebruikt om de resultaten voor de paper te genereren
d = length(startpoints_c); % nbr of parameters
min_v = zeros(d,1);
max_v = min_v;
n_v = min_v;
trans_m = zeros(d,2);
valid_c = cell(d,1);
idx_c = valid_c;

for i = 1:d
    min_v(i) = startpoints_c{i}(1);
    max_v(i) = startpoints_c{i}(2);
    n_v(i) = startpoints_c{i}(3);
    trans_m(i,:) = [min_v(i) max_v(i)-min_v(i)];
    valid_c{i} = linspace(0,1,nbr_valid(i));
    idx_c{i} = linspace(0,1,n_v(i));
end

valid_m = combvec( valid_c{:})'; % alle combinaties op een grid
aantal_valid = length(valid_m);
start_m = combvec( idx_c{:})';

aantal_start = length(start_m);
%% Initiele basis maken
A_conc = matrices{1};
B_conc = matrices{2};

n = size(A_conc,1);
nA_sum = size(A_conc,2)/n;
nB_sum = size(B_conc,2)/n;

omega_added = zeros(aantal_start,d);
timings_deriv = 0;
teller_afgel = 0;
timings_eigenv = 0;
teller_eigenv = 0;
%% Initialisation subspace (line 4)
for i = 1:aantal_start
    % Define matrices
    omega = trans_m(:,1) + diag(start_m(i,:))*trans_m(:,2);
    input_cell = num2cell( omega); 
    [alpha, beta, dalpha_c, dbeta_c] = functions( input_cell{:});
    A = A_conc*kron(alpha, speye(n)); 
    B = B_conc*kron(beta, speye(n));
    dA_c = cellfun(@(x) A_conc*kron(x, speye(n)), dalpha_c, 'UniformOutput', false);
    dB_c = cellfun(@(x) B_conc*kron(x, speye(n)), dbeta_c, 'UniformOutput', false);
    
    
    if i == 1
        n = size(A,1);
        basis_help = zeros(n,0);
        if strcmp( stan_gen, 'gen')
            [~, eig_B, flag] = eigs(B, 1, 'sa');
            if flag == 1
                eig_B = min( real( eig(B)));
            end
        else
            eig_B = 1;
            B = speye(n);
        end
    end
    %% calculating 1st eigenvalue
    start_eigenv = tic;
    [V, Lambda1] = calcul_eigv( A, B, sm_big, stan_gen, opts_eigs);
     
    V = V*diag(1./sqrt(diag(V'*B*V)));

    %% calculating 2nd eigenvalue
    [V2_sm, Lambda2] = eig( V(:,2:end)'*A*V(:,2:end), 'vector');    
    if strcmp( sm_big, 'sm')
        [Lambda2, idx] = min(Lambda2);
    else
        [Lambda2, idx] = max(Lambda2);
    end
    V2 = V(:,2:end)*V2_sm(:,idx);
    
    end_eigenv = toc(start_eigenv);   
    timings_eigenv = timings_eigenv + end_eigenv;
    teller_eigenv = teller_eigenv + 1;
    omega_added(i,:) = omega';
    %% derivative + adding
    if abs(Lambda2-Lambda1) > 10^(-8) && strcmp(method, 'noderiv') == 0
        if strcmp(method, 'system')
            start_afgel = tic;
            dX_m = dx_full(A,B,Lambda1, V(:,1), dA_c, dB_c);
            end_afgel = toc(start_afgel);
            timings_deriv = timings_deriv + end_afgel;
            teller_afgel = teller_afgel + 1;
        else
            [L,U, P, Q] = lu(A);
            start_afgel = tic;
            dX_m = deriv_eigv(A,B, dA_c, dB_c, Lambda1,V, {L,U, P, Q}, stan_gen, tol_gmres);
            end_afgel = toc(start_afgel);
            timings_deriv = timings_deriv + end_afgel;
            teller_afgel = teller_afgel + 1;
        end
        if nbr_vect == 1
            basis_help = [basis_help, V(:,1), dX_m];
        else
            basis_help = [basis_help, V(:,1), V2, dX_m];
        end
    else
        if nbr_vect == 1
            basis_help = [basis_help, V(:,1)];
        else
            basis_help = [basis_help, V(:,1), V2];
        end
    end
end

basis = orth( basis_help);
omega_done = [];
Abasis_conc = A_conc*kron( speye(nA_sum), basis);
if strcmp(stan_gen, 'gen')
    Bbasis_conc =  B_conc*kron( speye(nB_sum), basis);
end
n_basis = size(basis,2);

max_bovengrens = 1;

error_idx_v = 1:aantal_valid; % sorting error
isfulfilled = zeros(aantal_valid,1);
bovengrens_v = ones(aantal_valid,1);
lambda1_m = zeros(aantal_valid,d+1);
while (max_bovengrens > tol) 
    if size(basis, 2) > max_vect
        break
    end
    max_bovengrens = -1;
    aantal_valid - sum(isfulfilled) % nbr of points that need to be processed
    processed_points = 0;
    
    % line 8 in algorithm
    for i = 1:aantal_valid
        j = error_idx_v(i);
        if max_bovengrens > bovengrens_v(j)
            break;
        end
        if isfulfilled(j) == 1
            continue
        end
        processed_points = processed_points + 1;
        omega = trans_m(:,1) + diag(valid_m(j,:))*trans_m(:,2);
        input_cell = num2cell( omega); 
        [alpha, beta] = functions( input_cell{:});
        Abasis = Abasis_conc*kron(alpha, speye(n_basis));
        Ared = basis'*Abasis;
        if strcmp(stan_gen, 'gen')
            Bbasis = Bbasis_conc*kron(beta, speye(n_basis));
            Bred = basis'*Bbasis;
        end

        %% Solving small eigenvalue problem
        if strcmp( stan_gen, 'gen')
            if strcmp(sm_big,'sm')
                [hatX, hatLambda] = eig(Ared, Bred, 'vector');
                if nbr_vect > 1
                    [hatLambda,idx] = sort(real(hatLambda));
                    hatLambda2 = hatLambda(2);
                else
                    [hatLambda,idx] = min(real(hatLambda));
                end
                hatLambda_extr = hatLambda(1);
            else               
                [hatX, hatLambda] = eig(Ared, Bred, 'vector');
                if nbr_vect > 1
                    [hatLambda,idx] = sort(real(hatLambda),'descend');
                    hatLambda2 = hatLambda(2);
                else
                    [hatLambda, idx] = max(real(hatLambda));
                end
                hatLambda_extr = hatLambda(1);
            end               
        else
            if strcmp(sm_big,'sm')
                [hatX, hatLambda] = eig(Ared,'vector');
                if nbr_vect > 1
                    [hatLambda,idx] = sort(real(hatLambda));
                    hatLambda2 = hatLambda(2);
                else
                    [hatLambda, idx] = min(real(hatLambda));
                end
                hatLambda_extr = hatLambda(1);
            else
                [hatX, hatLambda] = eig(Ared,'vector');
                if nbr_vect > 1
                    [hatLambda,idx] = sort(real(hatLambda),'descend');
                    hatLambda2 = hatLambda(2);
                else
                    [hatLambda,idx] = max(real(hatLambda));
                end
                hatLambda_extr = hatLambda(1);
            end
        end
        lambda1_m(j,:) = [omega',hatLambda_extr];
        hatX_extr = hatX(:, idx(1));
        
        %% Determine residual and upper bound 
        % We assume that A and B depend in an affine way on w
        if strcmp(stan_gen, 'stan')
            residu = norm( Abasis*hatX_extr-hatLambda_extr*basis*hatX_extr);
        else
            residu = norm( Abasis*hatX_extr-hatLambda_extr*Bbasis*hatX_extr)/B_norm(basis*hatX_extr, B);
        end
        if nbr_vect == 1
            bovengrens_v(j) = min( [bovengrens_v(j),residu/sqrt(eig_B)]);
        else
            bovengrens_v(j) = min( [bovengrens_v(j),residu/sqrt(eig_B), residu^2/( eig_B*(hatLambda2-hatLambda_extr))]);
        end

        if bovengrens_v(j) < tol
            isfulfilled(j) = 1;
        else
            if bovengrens_v(j) > max_bovengrens
                max_bovengrens = bovengrens_v(j);
                idx_max = j;
            end
        end
    end
    processed_points
    [~, error_idx_v] = sort(bovengrens_v, 'descend');  
    if max_bovengrens > tol
        % Adding vectors to the subspace (line 15)
        omega = trans_m(:,1) + diag(valid_m(idx_max,:))*trans_m(:,2);
        input_cell = num2cell( omega); 
        [alpha, beta, dalpha_c, dbeta_c] = functions( input_cell{:});
        A = A_conc*kron(alpha, speye(n)); 
        B = B_conc*kron(beta, speye(n));
        dA_c = cellfun(@(x) A_conc*kron(x, speye(n)), dalpha_c, 'UniformOutput', false);
        dB_c = cellfun(@(x) B_conc*kron(x, speye(n)), dbeta_c, 'UniformOutput', false);
        
        
        start_eigenv = tic;
        [V, Lambda1] = calcul_eigv( A, B, sm_big, stan_gen, opts_eigs);
        V = V*diag(1./sqrt(diag(V'*B*V)));
        %% calculating 2nd eigenvalue
        [V2_sm, Lambda2] = eig( V(:,2:end)'*A*V(:,2:end), 'vector');    
        if strcmp( sm_big, 'sm')
            [Lambda2, idx] = min(Lambda2);
        else
            [Lambda2, idx] = max(Lambda2);
        end
        V2 = V(:,2:end)*V2_sm(:,idx);
        end_eigenv = toc(start_eigenv);
        timings_eigenv = timings_eigenv + end_eigenv;
        teller_eigenv = teller_eigenv + 1;
    
        %% Adding vectors to the subspace (line 19)
        if abs(Lambda2-Lambda1) > 10^(-8) && strcmp(method, 'noderiv') == 0
            if strcmp(method, 'system')
                start_afgel = tic;
                dX_m = dx_full(A,B,Lambda1, V(:,1), dA_c, dB_c);
                end_afgel = toc(start_afgel);
                timings_deriv = timings_deriv + end_afgel;
                teller_afgel = teller_afgel + 1;
            else
                [L,U, P, Q] = lu(A);
                start_afgel = tic;
                dX_m = deriv_eigv(A,B, dA_c, dB_c, Lambda1,V, {L,U,P,Q}, stan_gen, tol_gmres);
                end_afgel = toc(start_afgel);
                timings_deriv = timings_deriv + end_afgel;
                teller_afgel = teller_afgel + 1;
            end
            if nbr_vect == 1
                basis = [basis, V(:,1), dX_m]; % line 22 in algorithm
            else
                basis = [basis, V(:,1), V2, dX_m];
            end
        else
            if nbr_vect == 1
                basis = [basis, V(:,1)];
            else
                basis = [basis, V(:,1), V2];
            end
        end
        basis = orth(basis);
        % enkel voor sum
        Abasis_conc = A_conc*kron( speye(nA_sum), basis);
        if strcmp(stan_gen, 'gen')
            Bbasis_conc = B_conc*kron( speye(nB_sum), basis);
        end
        n_basis = size(basis,2);
        omega_added = [omega_added; omega'];
        lambda1_m(idx_max,:) = [omega', Lambda1];
        isfulfilled(idx_max) = 1;
        bovengrens_v(idx_max) = 0;
    end 
end
timings_deriv = timings_deriv/teller_afgel;
timings_eigenv = timings_eigenv/teller_eigenv;
if any(bovengrens_v > tol)
    display('Not all parameter values have a small enough error')
end