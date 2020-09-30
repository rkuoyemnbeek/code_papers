function [D1, D2] = D1_D2_maker( m, type)
%D1_D2_MAKER - This function makes the matrices needed for the definition of 
% the matrices in the convection-diffusion example. 
% Remark: this function maakt de matrices h*D1 en h*D2.

%
% INPUT:
%   (*) m = size van de matrices
%   (*) type = type boundary condities: 0: Dirichlet, 1: Neumann

%   Author: Koen Ruymbeek   
%   Address: Celestijnenlaan 200A, 3001 Leuven
%   email: koen.ruymbeek@cs.kuleuven.be
%   Website: https://www.kuleuven.be/wieiswie/nl/person/00114268
%   Date: 11-Oct-2019; Last revision: 11-Oct-2019
%
%   Copyright (c) 2019, Author
%   All rights reserved.

h = 1/(m+1);
D1 = 1/(2*h)*(spdiags( ones(m,1),1,m,m) - spdiags(ones(m,1),-1,m,m));
D2 = 1/h^2*(-2*spdiags( ones(m,1), 0,m,m) + spdiags( ones(m,1),1,m,m)+ spdiags( ones(m,1),-1,m,m)); 
if type == 1
    D1(end,1) = 1/(2*h);
    D1(1,end) = -1/(2*h);
    D2(end,1) = 1/h;
    D2(1,end) = 1/h;
end
