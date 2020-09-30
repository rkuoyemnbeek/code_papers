function dX_m = dx_full(A,B,Lambda1,x1, dA_c, dB_c)
 
d = length(dA_c);
n = size(A,1);
system_b = zeros(n+1,d);
system_A = [Lambda1*B-A B*x1; x1'*B 0];

for i = 1:d
    system_b(:,i) = [(dA_c{i}-Lambda1*dB_c{i})*x1; -(x1'*dB_c{i}*x1)/2];
end
solution = system_A\system_b;
dX_m = solution(1:end-1,:);
%dX_c = mat2cell(solution(1:end-1,:), n, ones(1,d));

