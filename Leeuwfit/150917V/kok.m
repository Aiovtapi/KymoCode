% Data

% A = [109 104 102
%     104 105 104
%     102 105 102
%     105 102 102
%     102 102 104
%     106 107 105];
% 
% B = [105 102
%     102 102];

% Engine

i=3;
j=1;

ROIr{i,j}=round(ROI{i,j});
ydatacrpdR1r{i,j}=round(ydatacrpdR1{i,j});

[iB,jB]=ndgrid(1:size(ROIr{i,j},1),1:size(ROIr{i,j},2));
is=sub2ind(size(ydatacrpdR1r{i,j}),iB,jB);
idx=bsxfun(@plus,(0:numel(ydatacrpdR1r{i,j})-is(end)),is(:));
imatch=find(all(bsxfun(@eq,ydatacrpdR1r{i,j}(idx),ROI{i,j}(:)),1));

% Display

for k=1:length(imatch)
    [i,j]=ind2sub(size(ydatacrpdR1r{i,j}),imatch(k));
    fprintf('upper-left corner match @ydatacrpd{i,j}(%d,%d)\n', i,j)
end
