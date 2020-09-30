clear all
close all
rng('default')
rng(0)

n1 = 30;
n2 = 30;
nu_v = linspace(0.9, 1.1, n1)';
xi_v = linspace(0.1, 0.2, n2)';
[W1, W2] = meshgrid( xi_v, nu_v);
    
for j = [1 3]
    [A0, A1,A2, f1, f2] = Delayed_Oscillator(500, j);
    A_c = {A0, A1, A2};
    n = size(A0,1);

% We want to approximate the rightmost eigenvalue for all parameter
% values (nu_v(i), xi_v(j)), i = 1,...,n1, j = 1,...,n2.
% We make tensor F_t with all function values, 
% so F_t(i0,i1,i2) = f_i0( nu_v(i1), xi_v(i2)) in this particular case
    F_t = zeros(3, n1, n2);
    F_t(1,:,:) = kron( ones(n1,1), ones(1,n2));
    F_t(2,:,:) = kron( f1(nu_v), ones(1,n2));
    F_t(3,:,:) = kron( -2*nu_v, xi_v');

    x0 = randn(n,1);
    max_rank = 50;
    max_iter = 100;

    [V, X_t, lambda, residu_v] = PSIRA_dD( A_c, F_t, x0, 10^(-3), 10^(-3),max_rank, max_iter, 10^(-10));
    figure(3)
    semilogy( residu_v/(sqrt(n1*n2)*normest(A0)))
    xlabel('iteration k')
    ylabel('\text{mean} \norm{\VM r(\omega_i)}/\norm{A_0}')
       
%  The approximate eigenvector of A(w^1_i1, w^2_i2, \hdots, w^d_id) is V*X_t(:,i1,...,id)
%  the resulting eigenvalues are contained in lambda.
%  Now the spectral abscissa for A(w) are (approx) the same as the spectral
%  abcissa for V'*A(w)*V. For new samples, it is enough to calculate the
%  spectral abscissa of V'*A(w)*V (which is a much smaller eigenvalue
%  problem.

    figure(1)
    surf( W1, W2, real( reshape( lambda, n1, n2)), 'EdgeColor', 'None')
    set(gcf,'position',[10,390,550/1.6,400/1.6])
    xlabel( '\omega_2')
    ylabel( '\omega_1')
    zlabel( '\Re(\lambda)')
    view(-63, 42)


end