function W_Z_tt = W_Z_maker_rank1( Z_tt, F_tt, U, q_m)
%W_Z_MAKER - Makes the TT-decomposition of W_Z_tt directly using the TT-decomposition 
% of Z_tt and F_tt. More details can be found in section 5.1 in the paper.
% This code is especially optimized for the case when all matrices are rank
% one except matrix A_0



%   Author: Koen Ruymbeek   
%   Address: Celestijnenlaan 200A, 3001 Leuven
%   email: koen.ruymbeek@cs.kuleuven.be
%   Website: https://www.kuleuven.be/wieiswie/nl/person/00114268
%   Date: 24-Sep-2019; Last revision: 7-febr-2020
%
%   Copyright (c) 2019, Author
%   All rights reserved.


d = Z_tt.d;
rz_v = Z_tt.r;
nz_v = Z_tt.n;
pz_v = Z_tt.ps;
coreZ = Z_tt.core;
rf_v = F_tt.r;
nf_v = F_tt.n;
m = nf_v(1)-1;
pf_v = F_tt.ps;
coreF = F_tt.core;
%% initialise X_Z_tt

rxz_v = rz_v.*rf_v;
nxz_v = [nz_v(1)+m; nz_v(2:end)];
coreXZ = zeros( [rxz_v(1:end-1).*nxz_v]'*rxz_v(2:end),1);
pxz_v = cumsum( [1; rxz_v(1:end-1).*nxz_v.*rxz_v(2:end)]);

%% for i from 1 to d-1
for i = 1:d-1
    GZ = reshape( coreZ(pz_v(i):pz_v(i+1)-1), rz_v(i)*nz_v(i), rz_v(i+1));
    GF = reshape( coreF(pf_v(i):pf_v(i+1)-1), rf_v(i)*nf_v(i), rf_v(i+1));
    GXZ = zeros( rxz_v(i)*nxz_v(i), rxz_v(i+1));
    if i == 1        
        GXZ(1:nz_v(1),:) = kron( GF(1,:), GZ);
        GXZ_help = q_m'*U*GZ;
        for j = 1:m
            GXZ(nz_v(1)+j,:) = kron( GF(j+1,:), GXZ_help(j,:));
        end
    else        
        for j = 1:nxz_v(i)
            GXZ( (j-1)*rxz_v(i)+1:j*rxz_v(i),:) = kron( GF( (j-1)*rf_v(i)+1:j*rf_v(i),:), GZ( (j-1)*rz_v(i)+1:j*rz_v(i),:));
        end
    end
    coreXZ(pxz_v(i):pxz_v(i+1)-1) = GXZ(:);
end

%% special case if i = d
GXZ = zeros( rxz_v(end-1), nxz_v(end));
GZ = reshape( coreZ(pz_v(end-1):pz_v(end)-1), rz_v(end-1), nz_v(end));
GF = reshape( coreF(pf_v(end-1):pf_v(end)-1), rf_v(end-1), nf_v(end));
for j = 1:rf_v(end-1)
    GXZ( (j-1)*rz_v(end-1)+1:j*rz_v(end-1),:) = GZ*diag( GF(j,:));
end
coreXZ(pxz_v(end-1):pxz_v(end)-1) = GXZ(:);
%% put all variables in the Tensor-Train Toolbox.
W_Z_tt = tt_tensor;
W_Z_tt.d = d;
W_Z_tt.r = rxz_v;
W_Z_tt.n = nxz_v;
W_Z_tt.core = coreXZ;
W_Z_tt.ps = pxz_v;
