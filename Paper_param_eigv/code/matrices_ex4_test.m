function varargout = matrices_ex4_test(x, y)
% matrices horende bij voorbeeld 4 zie https://arxiv.org/pdf/1504.06096.pdf
theta = [4/y; y/4; x*y/4];
if nargout > 2
    varargout = {theta, 1, {[0,0,y/4]', [-4/y^2, 1/4, x/4]'}, {0,0} };
else
    varargout = {theta, 1};
end