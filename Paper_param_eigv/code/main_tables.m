% main_tables: implementation of algorithm 1 of the paper. 
% We generate here the results of Table 1 and Table 2 

clear all
opts.issym = 1;
opts.spdB = 1;
opts1 = opts;
opts1.p = 70;
opts1.maxit = 500;

max_vect = 400;
for nbr_vb = 1 % nbr_vb1: table 1 in paper, nbr_vb2: table 2 in paper
    table = zeros(5,4);
    if nbr_vb == 1
        startpunten_c = {[-1/10, 1/10, 3], [2/10, 3/10, 3]};
        trainingset = [25 40];
        load('ex3_matrices.mat'); % 
        functions = @(x,y) matrices_ex3_test(x,y);
    elseif nbr_vb == 2
        startpunten_c = {[0.02 0.5 4], [2, 8, 4]};
        trainingset = [25 40];
        load('ex4_matrices.mat');
        functions = @(x,y) matrices_ex4_test(x,y);
    end

    for p = 1:2
        for k = 1:2
            A_conc = cell2mat(A_c);
            B_conc = cell2mat(B_c);
            matrices = {A_conc, B_conc};
            stan_gen_c = {'gen', 'gen'};
            if p == 1
                 start_duration = tic;
                [basis, omega_added,  Lambda_schat_valid_1, timings_afgel, timings_eigenv] = proj_algo_sym_RedBasis_satur(matrices, functions, max_vect, k, startpunten_c, trainingset, 10^(-5), 'sm', stan_gen_c{nbr_vb}, 'noderiv', 10^(-6), opts);            
                 duration = toc(start_duration);
            elseif p == 2
                 start_duration = tic;
                [basis, omega_added,  Lambda_schat_valid_2, timings_afgel, timings_eigenv] = proj_algo_sym_RedBasis_satur(matrices, functions, max_vect, k, startpunten_c, trainingset, 10^(-5), 'sm', stan_gen_c{nbr_vb}, 'system', 10^(-6), opts);            
                 duration = toc(start_duration);
            end
            d = length(startpunten_c);
            table(1,(p-1)*2+k) = size(basis,2);
            table(2,(p-1)*2+k) = length(omega_added);
            table(3,(p-1)*2+k) = duration;
            table(4,(p-1)*2+k) = timings_afgel/d;
            table(5,(p-1)*2+k) = timings_eigenv/k;
        end
    end
    table
end