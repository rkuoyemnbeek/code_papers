function [basis, omega_added,  Lambda_schat_valid, timings_afgel, timings_eigenv] = proj_algo_sym_RedBasis(matrices, functions, max_vect, aantal_vect, startpunten_c, n_valid, tol, sm_big, stan_gen, method, tol_gmres, opts_eigs)

d = length(startpunten_c); % aantal parameters
min_v = zeros(d,1);
max_v = min_v;
n_v = min_v;
trans_m = zeros(d,2);
valid_c = cell(d,1);
idx_c = valid_c;

for i = 1:d
    min_v(i) = startpunten_c{i}(1);
    max_v(i) = startpunten_c{i}(2);
    n_v(i) = startpunten_c{i}(3);
    trans_m(i,:) = [min_v(i) max_v(i)-min_v(i)];
    valid_c{i} = linspace(0,1,n_valid(i));
    idx_c{i} = linspace(0,1,n_v(i));
end

valid_m = combvec( valid_c{:})'; % alle combinaties op een grid
aantal_valid = length(valid_m);
Lambda_schat_valid = zeros( aantal_valid, d+1);
start_m = combvec( idx_c{:})';

aantal_start = length(start_m);
Lambda_schat_start = zeros( aantal_start, d+1);
%% Initiele basis maken
A_conc = matrices{1};
B_conc = matrices{2};

n = size(A_conc,1);
nA_sum = size(A_conc,2)/n;
nB_sum = size(B_conc,2)/n;

omega_added = zeros(aantal_start,d);
timings_afgel = 0;
teller_afgel = 0;
timings_eigenv = 0;
teller_eigenv = 0;
%% Initialisatie deelruimte 
for i = 1:aantal_start
    % Definieer matrices
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
    [V, Lambda1, R_B, perm_B] = calcul_eigv( A, B, sm_big, stan_gen, opts_eigs);
     
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
    
    cholB_c = {R_B, perm_B};
    Lambda_schat_start(i,:) = [omega' Lambda1];
    omega_added(i,:) = omega';
    %% derivative + adding
    if abs(Lambda2-Lambda1) > 10^(-8) && strcmp(method, 'noderiv') == 0
        if strcmp(method, 'system')
            start_afgel = tic;
            dX_m = dx_full(A,B,Lambda1, V(:,1), dA_c, dB_c);
            end_afgel = toc(start_afgel);
            timings_afgel = timings_afgel + end_afgel;
            teller_afgel = teller_afgel + 1;
        else
            [L,U, P, Q] = lu(A);
            start_afgel = tic;
            dX_m = deriv_eigv(A,B, dA_c, dB_c, Lambda1,V, {L,U, P, Q}, stan_gen, tol_gmres);
            end_afgel = toc(start_afgel);
            timings_afgel = timings_afgel + end_afgel;
            teller_afgel = teller_afgel + 1;
        end
        if aantal_vect == 1
            basis_help = [basis_help, V(:,1), dX_m];
        else
%             basis_help = [basis_help, V(:,1), dX_m];
            basis_help = [basis_help, V(:,1), V2, dX_m];
        end
    else
        if aantal_vect == 1
            basis_help = [basis_help, V(:,1)];
        else
%             basis_help = [basis_help, V(:,1)];
            basis_help = [basis_help, V(:,1), V2];
        end
    end
end
% basis = 0;
% omega_done = 0;
% timings_afgel/teller_afgel
% basis_no_orth = basis_help;
basis = orth( basis_help);
omega_done = [];
% als probleem kan geschreven worden als som.
Abasis_conc = A_conc*kron( speye(nA_sum), basis);
if strcmp(stan_gen, 'gen')
    Bbasis_conc =  B_conc*kron( speye(nB_sum), basis);
end
n_basis = size(basis,2);

max_error = 1;
teller = 1;
% load('lambda_ex3.mat')
% Lambda_m = Lambda_m(:);
while (max_error > tol) & (isempty(valid_m) == 0)
    if size(basis, 2) > max_vect
        break
    end
    aantal_valid = size(valid_m,1)
    bovengrens_v = zeros(aantal_valid,1);
    lambda1_m = zeros(aantal_valid,d+1);
%     start_valid = tic;
    for i = 1:aantal_valid
        %% bepalen van matrix
        omega = trans_m(:,1) + diag(valid_m(i,:))*trans_m(:,2);
        input_cell = num2cell( omega); 
        [alpha, beta] = functions( input_cell{:});
%         A = A_conc*kron(alpha, speye(n)); 
%         B = B_conc*kron(beta, speye(n));
%          Lambda1 = eigs(A,B,1, 'sm');
%         Lambda1 = lambda_v(i);
        Abasis = Abasis_conc*kron(alpha, speye(n_basis));
        Ared = basis'*Abasis;
        if strcmp(stan_gen, 'gen')
            Bbasis = Bbasis_conc*kron(beta, speye(n_basis));
            Bred = basis'*Bbasis;
        end

        %% klein eigenwaardeprobleem oplossen
        if strcmp( stan_gen, 'gen')
            if strcmp(sm_big,'sm')
                [hatX, hatLambda] = eig(Ared, Bred, 'vector');
                if aantal_vect > 1
                    [hatLambda,idx] = sort(real(hatLambda));
                    hatLambda2 = hatLambda(2);
                else
                    [hatLambda,idx] = min(real(hatLambda));
                end
                hatLambda_extr = hatLambda(1);
            else               
                [hatX, hatLambda] = eig(Ared, Bred, 'vector');
                if aantal_vect > 1
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
                if aantal_vect > 1
                    [hatLambda,idx] = sort(real(hatLambda));
                    hatLambda2 = hatLambda(2);
                else
                    [hatLambda, idx] = min(real(hatLambda));
                end
                hatLambda_extr = hatLambda(1);
            else
                [hatX, hatLambda] = eig(Ared,'vector');
                if aantal_vect > 1
                    [hatLambda,idx] = sort(real(hatLambda),'descend');
                    hatLambda2 = hatLambda(2);
                else
                    [hatLambda,idx] = max(real(hatLambda));
                end
                hatLambda_extr = hatLambda(1);
            end
        end
        lambda1_m(i,:) = [omega',hatLambda_extr];
        hatX_extr = hatX(:, idx(1));
        
        %% residu en bovengrens bepalen % we gaan uit dat A(w) = som(.) en B(w) = som(.)
        if strcmp(stan_gen, 'stan')
            residu = norm( Abasis*hatX_extr-hatLambda_extr*basis*hatX_extr);
        else
            residu = norm( Abasis*hatX_extr-hatLambda_extr*Bbasis*hatX_extr)/B_norm(basis*hatX_extr, B);
%            residu = norm( Abasis*hatX_extr/B_norm(basis*hatX_extr, B)-hatLambda_extr*Bbasis*hatX_extr/B_norm(basis*hatX_extr, B));

        end
        if aantal_vect == 1
            bovengrens_v(i) = residu/sqrt(eig_B);
        else
            bovengrens_v(i) = min( residu/sqrt(eig_B), residu^2/( eig_B*(hatLambda2-hatLambda_extr)));
        end         
    end
%     end_valid = toc(start_valid);
%     timings_valid = timings_valid + end_valid;
    idx_tol = find( bovengrens_v < tol);
    Lambda_schat_valid(teller:(teller+length(idx_tol)-1),:) = lambda1_m(idx_tol,:);
    teller = teller + length(idx_tol);
    bovengrens_v(idx_tol) = [];
    valid_m(idx_tol,:) = [];
    [max_error, idx_max] = max(bovengrens_v);
    if max_error > tol 
        % matrices bepalen
        omega = trans_m(:,1) + diag(valid_m(idx_max,:))*trans_m(:,2);
        input_cell = num2cell( omega); 
        [alpha, beta, dalpha_c, dbeta_c] = functions( input_cell{:});
        A = A_conc*kron(alpha, speye(n)); 
        B = B_conc*kron(beta, speye(n));
        dA_c = cellfun(@(x) A_conc*kron(x, speye(n)), dalpha_c, 'UniformOutput', false);
        dB_c = cellfun(@(x) B_conc*kron(x, speye(n)), dbeta_c, 'UniformOutput', false);
        
        
        start_eigenv = tic;
        [V, Lambda1, R_B, perm_B] = calcul_eigv( A, B, sm_big, stan_gen, opts_eigs);
        V = V*diag(1./sqrt(diag(V'*B*V))); % normaliseren;
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


        %% afgeleide bepalen + toevoegen vectoren aan basis
        if abs(Lambda2-Lambda1) > 10^(-8) && strcmp(method, 'noderiv') == 0
            if strcmp(method, 'system')
%                dX_m = dx_full(A,B,Lambda1, V(:,1), dA_c, dB_c);
                start_afgel = tic;
                dX_m = dx_full(A,B,Lambda1, V(:,1), dA_c, dB_c);
                end_afgel = toc(start_afgel);
                timings_afgel = timings_afgel + end_afgel;
                teller_afgel = teller_afgel + 1;
            else
%                dX_m = deriv_eigv(A,B, dA_c, dB_c, Lambda1,V, cholB_c, stan_gen, tol_gmres);
%                 [L,U] = lu(A);
                [L,U, P, Q] = lu(A);

                start_afgel = tic;
                dX_m = deriv_eigv(A,B, dA_c, dB_c, Lambda1,V, {L,U,P,Q}, stan_gen, tol_gmres);
                end_afgel = toc(start_afgel);
                timings_afgel = timings_afgel + end_afgel;
                teller_afgel = teller_afgel + 1;
            end
            if aantal_vect == 1
                basis = [basis, V(:,1), dX_m];
%                  basis_no_orth = [basis_no_orth, V(:,1), dX_m];
            else
                basis = [basis, V(:,1), V2, dX_m];
%                 basis_no_orth = [basis_no_orth, V(:,1), V2, dX_m];
            end
        else
            if aantal_vect == 1
                basis = [basis, V(:,1)];
%                 basis_no_orth = [basis_no_orth, V(:,1)];
            else
               basis = [basis, V(:,1), V2];
%                 basis_no_orth = [basis_no_orth, V(:,1), V2];
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
        valid_m(idx_max,:) = [];
        Lambda_schat_valid(teller,:) = [omega', Lambda1];
        teller = teller + 1;
    end 
end
timings_afgel = timings_afgel/teller_afgel;
timings_eigenv = timings_eigenv/teller_eigenv;