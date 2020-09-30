function [U]=round_rightleft_1mode(W_U, W_Z, erel)
% This is an adaptation from the function round in TT-toolbox. As we only
% need a basis for the columns space I adapt the code.

% subfunctions: lowrank_R

% Approximate TT-tensor with another one with specified accuracy
%   u=ROUND(TT,EPS) Approximate TT-tensor with relative accuracy EPS
%
%   u=ROUND(TT,EPS,RMAX) Approximate TT-tensor with relative accuracy 
%   EPS and maximal rank RMAX. RMAX can be array of ranks or a number
%
%
% TT-Toolbox 2.2, 2009-2012
%
%This is TT Toolbox, written by Ivan Oseledets et al.
%Institute of Numerical Mathematics, Moscow, Russia
%webpage: http://spring.inm.ras.ru/osel

d=W_Z.d;
n=W_Z.n;
r=W_Z.r;
pos=W_Z.ps;
cr=W_Z.core;
nrm=zeros(d,1);
core0=cr(pos(end-1):end);
%Orthogonalization from left-to-tight
for i=d:-1:2
    core0 = reshape(core0,[r(i),r(i+1)*n(i)]);
    [core0_help, ru_help] = qr( core0',0);
    core0 = core0_help';
    ru = ru_help'; 
    nrm(i-1)=norm(ru,'fro');
    
   if (nrm(i-1)~=0)
    ru=ru./nrm(i-1);
   end
   core1 = cr( pos(i)-r(i-1)*n(i-1)*r(i):pos(i)-1);
   core1=reshape(core1,[n(i-1)*r(i-1),r(i)]);
   core1=core1*ru;
   r(i) = size(core0,1);
   core0 = core1;
end
real_core0 = real(core0);
idx_imag = any(imag(core0));
core0 = [real_core0, imag(core0(:,idx_imag))];
U = lowrank_R( W_U, core0, erel);