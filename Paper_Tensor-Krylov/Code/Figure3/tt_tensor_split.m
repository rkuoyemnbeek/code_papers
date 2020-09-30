function [U0, tt] = tt_tensor_split( X, eps)
%TT_TENSOR_SPLIT - omzettting naar U0 en tt zodat U0 .1 tt = X
% X is full tensor, ipv deze eerst volledig om te zetten in
% TT-formaat en dan de basis voor de mode-1 tensor eruit te vissen, gaan we
% de code aanpassen zodat we rechtstreeks de U0 en tt hebben, zodat U0 .1 tt = X
% code verder aangepast zodat als X complex is, we in de 1ste iteratie 
% (om basis kolomruimte te bepalen)eerst X opsplitsen in reel en complex
% deel, dan svd bepalen en deze dan terug bij elkaar te zetten om dan
% verder te gaan

%
%   Other m-files required: TT-Toolbox-master

%
%   See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

%   Author: Koen Ruymbeek   
%   Address: Celestijnenlaan 200A, 3001 Leuven
%   email: koen.ruymbeek@cs.kuleuven.be
%   Website: https://www.kuleuven.be/wieiswie/nl/person/00114268
%   Date: 18-Apr-2019; Last revision: 18-Apr-2019
%
%   Copyright (c) 2019, Author
%   All rights reserved.

U0 = zeros(3);
tt = tt_tensor;
d = ndims( X);
n = size(X)';
r = ones(d+1,1);
core = zeros( prod(n),1);
c = X;
pos=1;
ep=eps/sqrt(d-1);

for i=1:d-1
    m=n(i)*r(i);
    c=reshape(c,[m,numel(c)/m]);

    if i == 1
        idx_imag = any( imag(c));
        c_schat_imag = imag(c(:,idx_imag));    
        c = [real(c), c_schat_imag];
        [U,Z] = lowrank_X( c, ep);
        r(i+1) = size(U,2);
        nbr_complex = length(find( idx_imag));
        complex_v_m = Z(:,size(Z,2)-nbr_complex+1:end);
        c = Z(:,1:size(Z,2)-nbr_complex);
        c(:,idx_imag) = c(:,idx_imag) + complex_v_m*1i;

        I = zeros( r(2)^2,1);
        I( 1:r(2)+1:r(2)^2) = 1;
        core(1:r(2)^2, 1) = I; % dit komt neer op I = eye(r(2)) en dan I(:)
        U0 = U;
        pos=pos+r(i+1)*r(i+1);
    else
        [u,s,v]=svd(c,'econ');
        s=diag(s); r1=my_chop2(s,ep*norm(s));
        u=u(:,1:r1); s=s(1:r1);
        r(i+1)=r1;    
        v=v(:,1:r1);
        v=v*diag(s); c=v';
        core(pos:pos+r(i)*n(i)*r(i+1)-1, 1)=u(:);
        pos=pos+r(i)*n(i)*r(i+1);
    end      
end
core(pos:pos+r(d)*n(d)*r(d+1)-1)=c(:);
core=core(1:pos+r(d)*n(d)*r(d+1)-1);
n = [r(2); n(2:end)];
ps=cumsum([1;n.*r(1:d).*r(2:d+1)]);
tt.d=d;
tt.n=n;
tt.r=r;
tt.ps=ps;
tt.core=core;
% tt.over=0; 

end
