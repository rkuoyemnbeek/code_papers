close all
clear all

path_image_folder = '~/';
file_name = mfilename('fullpath');
date_name = date;
tikzcomment = sprintf('This file was created in matlab file %s on %s', file_name, date_name);

% Make sure you add the TT-Toolbox to your path.
rng( 'default');

%% initialize all matrices and vectors
c1 = 5;
c2 = 5;
d11 = -1;
d22 = -1;
d12 = 0;

m1 = 100;
m2 = 100;
n = m1*m2;
h1 = 1/(m1+1);
h2 = 1/(m2+1);

[D_1m1, D_2m1] = D1_D2_maker(m1,0);
[D_1m2, D_2m2] = D1_D2_maker(m2,0);

B1 = kron( D_1m2, D_1m1);
f1 = @(c1, c2, d11, d22, d12) d12;
B2 = kron( eye(m2), D_1m1);
f2 = @(c1, c2, d11, d22, d12) c1;
B3 = kron( eye(m2), D_2m1);
f3 = @(c1, c2, d11, d22, d12) d11;
B4 = kron( D_1m2, eye(m1));
f4 = @(c1, c2, d11, d22, d12) c2;
B5 = kron( D_2m2, eye(m1));
f5 = @(c1, c2, d11, d22, d12) d22;

A = f1(c1, c2, d11, d22, d12)*B1 + f2(c1, c2, d11, d22, d12)*B2 + f3(c1, c2, d11, d22, d12)*B3 + ...
    f4(c1, c2, d11, d22, d12)*B4 + f5(c1, c2, d11, d22, d12)*B5;

A1 = B2;
A2 = B4;
A3 = B3+B5;
grootte = [10, 30];
l = length(grootte);
tijd_2D = zeros(l,1);
tijd_1D = zeros(l,1);


tol_W = 10^(-3);
tol_X = 10^(-3);

A_c = {A1, A2, A3};
x0 = randn(n,1);
V = x0/norm(x0);  

if isempty(gcp('nocreate')) == 1 % check if there is a parallel pool, if not, we make one
    parpool( 'local') 
    maxNumCompThreads(5)
end

%% We calculate the eigenvalues over the parameter when n1=n2=n3=10 en n1=n2=n3=30

for i = 1:l
    c1_v = linspace(4,6, grootte(i))';
    c2_v = linspace(4,6, grootte(i))';
    d11_v = linspace(-1.1, -0.9, grootte(i))';
    w_c = {c1_v, c2_v, d11_v};
    n1 = length(c1_v);
    n2 = length(c2_v);
    n3 = length(d11_v);
    n_v = [n1,n2, n3];
    
    F_t = zeros(3,n1,n2,n3);
    F_t(1,:,:,:) = reshape(w_c{1}*ones(1,n2*n3), n1,n2,n3);
    G_thelp = w_c{2}*ones(1,n1*n3);
    F_t(2,:,:,:) = permute( reshape( G_thelp, n2,n1,n3),[2,1,3]);
    G_thelp = w_c{3}*ones(1,n1*n2);
    F_t(3,:,:,:) = permute( reshape( G_thelp, n3,n1,n2),[2,3,1]);
    
    F_tt = tt_tensor(F_t, 10^(-15));
    f1_v = F_t(1,:);
    f2_v = F_t(2,:);
    f3_v = F_t(3,:);
    max_rank = 100;
    
    eps = 10^(-9);
    [V, X_schat_t, lambda_schat_t, residu_v] = Arnoldi_dD( A_c, F_t, F_tt, w_c, V, tol_W, tol_X, max_rank,eps);
    save( strcat('conv_diff_3D_', num2str(i), '.mat'),'A_c','F_t','w_c', 'tol_W', 'tol_X', 'V','X_schat_t', 'lambda_schat_t', 'residu_v')
end

for i = 1:2
    load( strcat('conv_diff_3D_', num2str(i), '.mat'))
    figure(i)
    semilogy( residu_v/sqrt(grootte(i)^3))
    xlabel('iteration $k$')
    ylabel('$ \norm{\VM R}_F/\sqrt{n_1\hdots n_d}$')
end

