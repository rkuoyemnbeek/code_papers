function W_Z_tt = W_Z_maker( Z_tt, F_tt)
%X_Z_MAKER - Based on equation (19) in the paper pg 12, we make the
%TT-decomposition of W_Z. De tensors F_tt and Z_tt are de tensor F resp. Z
%in the paper.

%
%   Other m-files required: none
%   Subfunctions: none
%   MAT-files required: none
%

%   Author: Koen Ruymbeek   
%   Address: Celestijnenlaan 200A, 3001 Leuven
%   email: koen.ruymbeek@cs.kuleuven.be
%   Website: https://www.kuleuven.be/wieiswie/nl/person/00114268
%   Date: 24-Sep-2019; Last revision: 24-Sep-2019
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
pf_v = F_tt.ps;
coreF = F_tt.core;

%% initialiseren X_Z_tt
rxz_v = rz_v.*rf_v;
nxz_v = [nz_v(1)*nf_v(1); nz_v(2:end)];
coreXZ = zeros( [rxz_v(1:end-1).*nxz_v]'*rxz_v(2:end),1);
pxz_v = cumsum( [1; rxz_v(1:end-1).*nxz_v.*rxz_v(2:end)]);

%% voor i is 1 tot d-1
for i = 1:d-1
    GZ = reshape( coreZ(pz_v(i):pz_v(i+1)-1), rz_v(i)*nz_v(i), rz_v(i+1));
    GF = reshape( coreF(pf_v(i):pf_v(i+1)-1), rf_v(i)*nf_v(i), rf_v(i+1));
    if i == 1
        GXZ = kron( GF, GZ);
    else
        GXZ = zeros( rxz_v(i)*nxz_v(i), rxz_v(i+1));
        for j = 1:nxz_v(i)
            GXZ( (j-1)*rxz_v(i)+1:j*rxz_v(i),:) = kron( GF( (j-1)*rf_v(i)+1:j*rf_v(i),:), GZ( (j-1)*rz_v(i)+1:j*rz_v(i),:));
        end
    end
    coreXZ(pxz_v(i):pxz_v(i+1)-1) = GXZ(:);
end

%% voor i = d
GXZ = zeros( rxz_v(end-1), nxz_v(end));
GZ = reshape( coreZ(pz_v(end-1):pz_v(end)-1), rz_v(end-1), nz_v(end));
GF = reshape( coreF(pf_v(end-1):pf_v(end)-1), rf_v(end-1), nf_v(end));
for j = 1:rf_v(end-1)
    GXZ( (j-1)*rz_v(end-1)+1:j*rz_v(end-1),:) = GZ*diag( GF(j,:));
end
coreXZ(pxz_v(end-1):pxz_v(end)-1) = GXZ(:);

%% in TT-formaat steken
W_Z_tt = tt_tensor;
W_Z_tt.d = d;
W_Z_tt.r = rxz_v;
W_Z_tt.n = nxz_v;
W_Z_tt.core = coreXZ;
W_Z_tt.ps = pxz_v;
