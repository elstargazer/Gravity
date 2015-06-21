function a=VWVerticalRecursion(n, m, z, rsquared,RefRad,PreviousVW,PrePreviousVW)


a=(2*n-1).*z.*RefRad.*PreviousVW./(n-m)./rsquared-(n+m-1).*RefRad.*RefRad.*PrePreviousVW./rsquared./(n-m);