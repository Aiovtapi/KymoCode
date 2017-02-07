function [M,MSTD] = LionRowMS(K)
%LIONROWMEAN takes the nonzero mean of rows, but equates to zero if row
%only consists of zero.
% K is a cell with multiple elements (depending on spot number)
% M is the mean cell containing mean values per spot
% m indicates the column in M where the mean values should go
% m=1 intensity m=2 x pos m=3 x std m=4 y pos m=5 y std m=6 full cell
% intensity

Nspots=size(K,1);
BacLife=size(K{1},1);

for i=1:Nspots
    for j=1:BacLife
    if isempty(nonzeros(K{i}(j,:)))
        M{i}(j,1)=0;
        MSTD{i}(j,1)=0;
    else
        M{i}(j,1)=nanmean(nonzeros(K{i}(j,:)));
        MSTD{i}(j,1)=nanstd(nonzeros(K{i}(j,:)));
    end
    end
end

end

