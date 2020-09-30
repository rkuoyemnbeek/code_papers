function [A0, A1,A2, f1, f2] = Delayed_Oscillator(N, RE)
% I make from delay problem, matrices A0, A1 and A2 such that
% A(nu, xi) = A0 + f1(nu, xi) A1 + f2(nu, xi) A2 
n=2;    % Dimension of the system
h=2;    % Number of delays (you should count also the 0 delay)
D=2;    % Number of parameter which may vary
% N = 500;
% RE = 2;

% Different controllers - look at the definition of RE [1, Table 2]
if RE==1                    % SAE - SMOOTH CASE
    K=[0.2, 0.2];   
elseif RE==2                % MSSAE - Obtained by probabilistic stabilization kappa=0
    K=[0.5105, -9.1810e-2];
elseif RE==3                % MNSSAE  - Obtained by deterministic stabilization
    K=[0.6179, -7.1644e-3];
elseif RE==4                % Obtained by probabilistic stabilization kappa=10
    K=[0.45043,-0.1631];
elseif RE==5                % Obtained by probabilistic stabilization kappa=100
    K=[0.36328,-0.39879];    
elseif RE==6                % Obtained by worst-case stabilization 
    K=[0.44437,-0.15702];  
else
%    error('RE must be defined as the numbers 1,2 and 3 for SAE, MSSAE, MNSSAE respectively')
    K = [0 0];
end

% Definition of the Delay Differential Equation associated to the
% oscillator with feedback delay 
TAU=[0,1];      % Delays
% E=eye(n);        % Leading matrix NON SINGULAR
A=cell(1,h);
% A{1}= [   0,        1;...
%       -nu^2, -2*nu*xi];  
A{h}= [   0,    0;...
       K(1), K(2)]; 
   
   
x = cos(pi*(0:N)/N)';
x=(TAU(h)/2)*(x-1);     % SHIFT for the delay
c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
X = repmat(x,1,N+1);
dX = X-X';
D = (c*(1./c)')./(dX+(eye(N+1))); % off-diagonal entries
D = D - diag(sum(D')); % diagonal entries
DN=kron(D,eye(n));

A0=zeros(n*(N+1));
A0(1,n) = 1;
% See [3, example 5.1]
% AN(1:n,1:n)=A{1};
A0(1:n,end-n+1:end)=A{2};

A0(n+1:end,:)=DN(n+1:end,:);
% Definition of the eigenvalue problem
% A=A1;
A1=zeros(n*(N+1));
f1 = @(nu, xi) -nu.^2;
A1(1:n,1:n)=[0,0;1, 0];

f2 = @(nu, xi) -2*nu.*xi;
A2=zeros(n*(N+1));
A2(1:n,1:n)= [0,0;0, 1]; 



