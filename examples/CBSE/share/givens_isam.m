function [c, s] = givens_isam(x, y)
%Supplies givens matrix for a QR update step
% (slightly modified / optimized version of givens)
absx = abs(x);
if absx == 0.0
    c = 0.0; s = 1.0;
else
    nrm = norm([x y]);
    c = absx/nrm;
    s = x/absx*(conj(y)/nrm);
end