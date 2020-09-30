function varargout = matrices_ex3_test(x, y)
% matrices horende bij voorbeeld 3 zie https://arxiv.org/pdf/1504.06096.pdf
%% alpha
alpha = zeros(16,1);
alpha(1,1) = ((2*x + 1)*(y - 1))/((2*y - 1)*(y + 1));
alpha(2,1) = ((8*x^2)/3 + 1/2)/(y + 1) - (4*x^2)/(3*(2*y - 1));
alpha(3,1) = (2*x*(y - 1))/((2*y - 1)*(y + 1));
alpha(4,1) = (2*x^2 + 2/3)/(y + 1) - 1/(6*y - 3);
alpha(5,1) = (x + 1/2)/(y + 1);
alpha(6,1) = (y - 1)/((2*x + 1)*(2*y - 1)*(y + 1));
alpha(7,1) = - 1/(3*(y + 1)) - 1/(6*y - 3);
alpha(8,1) = -((2*x - 1)*(y - 1))/((2*y - 1)*(y + 1));
alpha(9,1) = -(x - 1/2)/(y + 1);
alpha(10,1) = -(y - 1)/((2*x - 1)*(2*y - 1)*(y + 1));
alpha(11,1) = -1/(2*(2*x - 1)*(y + 1));
alpha(12,1) = -x/((2*y - 1)*(y + 1));
alpha(13,1) = 1/(2*(2*x + 1)*(y + 1));
alpha(14,1) = 1/(2*(y + 1));
alpha(15,1) = 2/(3*(y + 1)) - 1/(3*(2*y - 1));
alpha(16,1) = x/(y + 1);

%% dalpha
% dx 
if nargout > 2
    dalpha_dx = zeros(16,1);

    dalpha_dx(1,1) = (2*(y - 1))/(2*y^2 + y - 1);
    dalpha_dx(2,1) = (8*x*(y - 1))/(2*y^2 + y - 1);
    dalpha_dx(3,1) =  (2*(y - 1))/((y + 1)*(2*y - 1));
    dalpha_dx(4,1) = (4*x)/(y + 1);
    dalpha_dx(5,1) = 1/(y + 1);
    dalpha_dx(6,1) = -(2*(y - 1))/((2*x + 1)^2*(2*y^2 + y - 1));
    dalpha_dx(7,1) = 0;
    dalpha_dx(8,1) =  -(2*(y - 1))/(2*y^2 + y - 1);    
    dalpha_dx(9,1) = -1/(y + 1);
    dalpha_dx(10,1) = (2*(y - 1))/((2*x - 1)^2*(2*y^2 + y - 1));    
    dalpha_dx(11,1) = 1/((2*x - 1)^2*(y + 1));
    dalpha_dx(12,1) = -1/((y + 1)*(2*y - 1));
    dalpha_dx(13,1) = -1/((2*x + 1)^2*(y + 1));
    dalpha_dx(14,1) = 0;
    dalpha_dx(15,1) = 0;
    dalpha_dx(16,1) = 1/(y + 1);

    % dy
    dalpha_dy = zeros(16,1);
    dalpha_dy(1,1) = -(2*(2*x + 1)*(y - 2)*y)/(2*y^2 + y - 1)^2;
    dalpha_dy(2,1) = (8*x^2)/(3*(2*y - 1)^2) - ((8*x^2)/3 + 1/2)/(y + 1)^2;
    dalpha_dy(3,1) =  -(4*x*(y - 2)*y)/(2*y^2 + y - 1)^2;
    dalpha_dy(4,1) = 6/(6*y - 3)^2 - (2*x^2 + 2/3)/(y + 1)^2;
    dalpha_dy(5,1) = -(x + 1/2)/(y + 1)^2;
    dalpha_dy(6,1) =  -(2*(y - 2)*y)/((2*x + 1)*(2*y^2 + y - 1)^2);    
    dalpha_dy(7,1) = (2*y^2 + 1)/(2*y^2 + y - 1)^2;   
    dalpha_dy(8,1) =  (2*(2*x - 1)*(y - 2)*y)/(2*y^2 + y - 1)^2;   
    dalpha_dy(9,1) = (2*x - 1)/(2*(y + 1)^2);    
    dalpha_dy(10,1) =  (2*(y - 2)*y)/((2*x - 1)*(2*y^2 + y - 1)^2);   
    dalpha_dy(11,1) = 1/(2*(2*x - 1)*(y + 1)^2);    
    dalpha_dy(12,1) = (4*x*y + x)/(2*y^2 + y - 1)^2;
    dalpha_dy(13,1) = -1/(2*(2*x + 1)*(y + 1)^2);

    dalpha_dy(14,1) = -1/(2*(y + 1)^2);
    dalpha_dy(15,1) = -(2*(y - 2)*y)/(2*y^2 + y - 1)^2;
    dalpha_dy(16,1) = -x/(y + 1)^2;   
end

if nargout < 3 
    varargout = {alpha, 1};
else
    varargout = {alpha, 1, {dalpha_dx, dalpha_dy}, {0, 0}};
end