function lmcosi_out = WeightSumExpansion(w,lmcosi_cell)

% all lmcosi's should have the same length;
% sum of weights should be 1

lmcosi_out = lmcosi_cell{1};
lmcosi_out(:,3:4)=0;

for i=1:numel(lmcosi_cell)
    lmcosi_out(:,3:4)=lmcosi_out(:,3:4) + w(i)*lmcosi_cell{i}(:,3:4);    
end