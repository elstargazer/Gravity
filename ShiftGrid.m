function ris=ShiftGrid(ri,lambdai,row,delta)


ris=interp1(lambdai(row,:),ri(row,:),lambdai(row,:)+delta,'linear');