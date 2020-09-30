function [A0, A1,A2, f1, f2] = Heating_Luca(N, RE)
% This function makes the example in Main_Figure5
% 
% References:
% [1] T. Vyhlidal, P. Zitek and K. Paulu, "Design, modelling and control
%     of the experimental heat transfer set-up", in Loiseau et al. (Eds.),
%     Topics in Time Delay Systems Analysis, Algorithms, and Control, Lecture
%     Notes in Control and Information Sciences, Springer, New York, 2009, pp.
%     303-314.
% [2] W. Michiels, T. Vyhlidal, P. Zitek, "Control design for time-delay
%     systems based on quasi direct pole placement", Journal of Process Control
%     20: 337-343 (2010).
% [3] L. Fenzi and W. Michiels, "Robust stability optimization for linear 
%     delay systems in a probabilistic framework", Linear Algebra Appl. 
%     526: 1-26, 2017.

% [4] L. Fenzi, "Looking for stability, advances on spectrum-based
%     stability and stabilization of uncertain linear time delay systems",
%     PhD thesis, KU Leuven, 2020.

n=10;    % Dimension of the system
h=8;    % Number of delays (you should count also the 0 delay)
D=2;    % Number of parameter which may vary

%% CONTROLLER
if RE==1        % DETERMINISTIC STABILIZATION
    K=[-2.0000e+00 , -9.8150e-01  , 1.0168e+00  , 3.6514e-01 , -9.2268e+00 ,  4.1945e+00 , -5.8979e+00  , 7.3821e-02 ,  9.9198e+00 ,  6.4099e-01];
elseif RE==2    % PROBABILISTIC STABILIZATION with 1 UNCERTAINTY kappa=0
    K=[-6.8858e-01, -1.2222e+00, -1.4796e+00,  2.1180e+00,  1.2478e+00, -8.1282e-02, -1.7132e+00, -1.2819e+00,-5.9722e-01,3.5813e-01];
elseif RE==3    % PROBABILISTIC STABILIZATION with 1 UNCERTAINTY kappa=1000
    K=[-6.8858e-01,  -1.2222e+00,  -1.4796e+00,   2.1180e+00,   1.2478e+00,  -8.1282e-02,  -1.7132e+00,  -1.2819e+00,  -5.9722e-01,   3.5813e-01];
elseif RE==4    % PROBABILISTIC STABILIZATION with 2 UNCERTAINTY kappa=0
    K=[-1.5001e+00, 5.9506e-01,-9.0441e-02,1.2190e-01,  -6.6453e-01,  -1.5547e+00,  -2.4243e-02,  -3.1575e-02,   1.2449e+00,   7.6500e-01];
elseif RE==5    % PROBABILISTIC STABILIZATION with 2 UNCERTAINTY kappa=1000
    K=[  -1.4997e+00,   5.9591e-01 , -8.9697e-02,   1.2136e-01,  -6.6513e-01,  -1.5549e+00,  -2.6209e-02,  -3.2337e-02,   1.2440e+00,   7.6385e-01];
end

%% Define Model Parameters
Trd1    = 5;
Krd1    = 0.973;
tau_rd1 = 15;

Trh         = 25;
Krh1        = 0.96;
Krh2        = 0.536;
tau_rh1     = 23;
tau_rh2     = 11.5;

Trd2        = 5;
Krd2        = 0.983;
tau_rd2     = 5;

Tre         = 3;
kappa_re    = 0.85;

Trd3        = 5;
Krd3        = 0.975;
tau_rd3     = 3;

Krc1        = 0.9;    
Krc2        = 0.043;   
tau_rc1     = 5;
tau_rc2     = 6;

Tll         = 63;
tau_ll1     = 29;
tau_ll2     = 7;

Tle         = 3;
kappa_le    = 0.85;

Tld         = 5;
Kld         = 0.97;
tau_ld      = 5;

Klc1        = 0.92;     
Klc2        = 0.042;    
tau_lc1     = 5;
tau_lc2     = 6;

%% UNCERTAIN PARAMETERS
% Temperature left cooler   (Tlc = 15+/- 10%)
% Temperature right cooler  (Trc = 17+/- 10%)

%% DISCRETE DELAY TERMS
tau=[0 3 5 6 7 15 23 29]; 
% LEADING MATRIX: E=eye(n);


%% FIRST SYSTEM MATRIX - A0
A0 = zeros(9,10); 
A0(1,1) = -1/Trd1; A0(2,2) = -1/Trh; A0(3,3) = -1/Trd2; A0(4,3) = (1-0.5*kappa_re)/Tre;
A0(4,4) = -(1+0.5*kappa_re)/Tre; A0(4,7) = kappa_re/2/Tre; A0(4,8) = kappa_re/2/Tre;
A0(5,5) = -1/Trd3; 
% A0(6,6) = -1/Trc;         % Uncertainty 2
A0(8,3) = kappa_le/2/Tle; A0(8,4) = kappa_le/2/Tle;
A0(8,7) = (1-0.5*kappa_le)/Tle; A0(8,8) = -(1+0.5*kappa_le)/Tle; A0(9,9) = -1/Tld;
% A0(10,10) = -1/Tlc;       % Uncertainty 1
A=[A0(1:5,:); 
   zeros(1,10);
   A0(7:9,:);
   zeros(1,10)];

dA=cell(1,2);
dA{1}=zeros(10); dA{1}(10,10)=-1;
dA{2}=zeros(10); dA{2}(6,6)=-1;
dfA=cell(1,2);
dfA{1}=@(omega1) 1/omega1;
dfA{2}=@(omega2) 1/omega2;

p=h-1; 
B=cell(p,1);    % Matrices at time tau>0
dB=cell(p,2);   % perturbation
dfB=cell(p,2);  % uncertain dependence

%% 2nd Matrix - A1
A1=zeros(10); A1(5,4) = Krd3/Trd3;
B{1}=A1;
dB{1,1}=zeros(10);
dB{1,2}=zeros(10);

%% 3rd Matrix - A2
A2 = zeros(9,10); A2(3,2) = Krd2/Trd2;
% A2(6,5) = Krc1/Trc; % Uncertainty 2
A2(9,8) = Kld/Tld;
% A2(10,9) = Klc1/Tlc; % Uncertainty 1
B{2}=[A2(1:5,:); 
    zeros(1,10);
    A2(7:9,:);
    zeros(1,10)];
dB{2,1}=zeros(10); dB{2,1}(10,9)= Klc1;
dB{2,2}=zeros(10); dB{2,2}(6,5) = Krc1;

%% 4th Matrix  - B1
%B1=zeros(10,1); B1(10)=Klc2/Tlc;
B{3}=zeros(10);
dB{3,1}=zeros(10); dB{3,1}(10,:)= Klc2*K;
dB{3,2}=zeros(10); 

%% 5th Matrix  - B2
%B2=zeros(10,1); B2(7)=1/Tll;
B{4}=zeros(10); B{3}(7,:)=K/Tll;
dB{4,1}=zeros(10); 
dB{4,2}=zeros(10); 

%% 6th Matrix  - A3
A3=zeros(10); A3(1,6) = Krd1/Trd1;
B{5}=A3;
dB{5,1}=zeros(10); 
dB{5,2}=zeros(10); 

%% 7th Matrix  - A4             DELAY 23
A4=zeros(10); A4(2,1) = Krh1/Trh;
B{6}=A4;
dB{6,1}=zeros(10); 
dB{6,2}=zeros(10); 

%% 8th Matrix  - A5             DELAY 29
A5=zeros(10); A5(7,7) = -1/Tll;
B{7}=A5;
dB{7,1}=zeros(10); 
dB{7,2}=zeros(10);

%% PARAMETER DEPENDENCIES
f1=@(omega1) 1/omega1;
f2=@(omega2) 1/omega2;

Id=eye(n); % Identity for the system

%% PIECEWISE MESH DEFINITION
dtau=diff(tau)/2;           %half-length of delay intervals
Mk=em1(2*dtau,p,N);         %piecewise spectral discretization of delay
                            %intervals: Mk(k)+1 nodes in
                            %[-tau_{k},-tau_{k-1}]
                            %em1: M points in total
                            %em2: M points on the largest interval
MMk=[0,cumsum(Mk)];         %cumulative Mk's
M=MMk(end);                 %M+1=total number of nodes in [-max(tau),0]


%% MATRIX BLOCK-DIAGONAL: DIFFERENTIATION PART
AM=zeros(M+1);                                %initialization
for k=1:p                                     %along [-tau_{k},-tau_{k-1}]
    OmegaMk=dtau(k)*cos((0:Mk(k))*pi/Mk(k))-...
        (tau(k)+tau(k+1))/2;                  %Chebyshev II nodes
    Dk=difmat(OmegaMk);                       %differentiation matrix
    AM(MMk(k)+2:MMk(k+1)+1,...
        MMk(k)+1:MMk(k+1)+1)=Dk(2:Mk(k)+1,:); %update
end                                           %
AM=kron(AM,Id);                               %block structure for systems

%% MATRIX 1ST BLOCK-ROW: SPLICING CONDITION
AM(1:n,1:n)=A;                    %update with current time matrix
for k=1:p                                     %along [-tau_{k},-tau_{k-1}]
    AM(1:n,MMk(k+1)*n+1:(MMk(k+1)+1)*n)=B{k}; %update
end                                           

%% OUTPUT
A0=AM;
A1=zeros(size(AM)); A2=zeros(size(AM));
A1(1:n,1:n)=dA{1};   A2(1:n,1:n)=dA{2};
for k=1:p                                     %along [-tau_{k},-tau_{k-1}]
    A1(1:n,MMk(k+1)*n+1:(MMk(k+1)+1)*n)=dB{k,1}; %update
    A2(1:n,MMk(k+1)*n+1:(MMk(k+1)+1)*n)=dB{k,2}; %update
end    


end

%% ADDITIONAL FUNCTIONS


function Mk=em1(d,p,M)
%EM1 piecewise spectral distribution.
%  Mk=EM1(d,p,M) distributes a total of sum(Mk) points over a grid of p
%  intervals of length d(1),...,d(p) proportionally to spectral accuracy,
%  i.e., (d(i)/Mk(i))^Mk(i) is almost constant for all i=1,...,p, under the
%  constraint that Mk(1),...,Mk(p) are the minimum integers s.t.
%  sum(Mk)>=M.
%
%  Version 1.0, may 28, 2014.

if isempty(d)
    Mk=0;
    return
end
Mk=zeros(1,p);
[~,ind]=max(d);
f=@(x,t,m)(t/x)^x-(d(ind)/m)^m;
m=ceil(M/p);
while sum(Mk)<M
    for k=1:p
        mm=m;
        while (f(mm,d(k),m)<0)&&(mm>=1)
            mm=mm-1;
        end
        Mk(k)=mm;
    end
    m=m+1;
end

end

function D=difmat(x)
%DIFMAT differentiation matrix on Chebyshev II nodes.
%  D=DIFMAT(x) returns the differentiation matrix D relevant to the N+1
%  Chebyshev II nodes x according to [3].
%
%  Original version from [3].

N=length(x)-1;
if N==0
    D=0;
    return
end
c=[2;ones(N-1,1);2].*(-1).^(0:N)';
X=repmat(x',1,N+1);
dX=X-X';
D=(c*(1./c)')./(dX+(eye(N+1)));
D=D-diag(sum(D'));
end