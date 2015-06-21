function [Vnew,Wnew]=VWDiagonalRecursion(m, x, y, rsquared, RefRad,V, W)

Vnew=(2*m-1).*RefRad.*(x.*V-y.*W)./rsquared;
Wnew=(2*m-1).*RefRad.*(x.*W+y.*V)./rsquared;