function [BGCell] = LionLoadBackground(utrack)
% This function translates the utrack output data to create a single cell
% with trajectory information in rows, instead of columns. This is the
% format:
%
% For each track:
% 
% create a cell containing matrices [X Y Z AMP Xstd Ystd Zstd Ampstd]

Nframes=size(utrack.Detection.localMaxima,1);
BGCell=cell(Nframes,1);

for i=1:Nframes
    
    Ndetections=size(utrack.Detection.localMaxima(i).cands,1);
    
    for j=1:Ndetections;
        
    BGCell{i}(j,1)=utrack.Detection.localMaxima(i).cands(j).IBkg*(2^16-1);
    
    end
end

% BGCell=BGCell(~cellfun('isempty',BGCell));

end