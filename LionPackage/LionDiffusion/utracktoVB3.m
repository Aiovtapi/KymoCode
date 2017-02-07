function [F] = utracktoVB3(tracksFinal,lengthlowerbound,lengthhigherbound,dim)
%UTRACKTOVB3 formats utrack trajectories to VB3 input format
% input format VB3 = [Amp, X, Y, Z], being column vectors.
% utrack: [x1 y1 z1 Amp1 xstd1 ystd1 zstd1 ampstd1 x2 y2 ....]

Ndatasets=size(tracksFinal,2);

for k=1:Ndatasets
    
    T=tracksFinal{k};
    Tracks=T.tracksFinal;
    
Ntracks=size(Tracks,1);
FinalTraj.Trajectories=cell(1,Ntracks);
umperpx=0.1589;
FinalTraj.dim=dim;

%extract from utrack file, to matrix with timesteps in rows and not being
%concatenated

for i=1:Ntracks
    TrackLength(i,1)=size(Tracks(i).tracksFeatIndxCG,2);
    for j=1:TrackLength(i,1)
    FinalTraj.Trajectories{i}(j,:)=Tracks(i).tracksCoordAmpCG(1+8*(j-1):8*j);
    end
    
    % adjust for dimensionality of the analysis.
    % could be determine from u-track's post-analysis, as it gives 0 for
    % dimensions not encounted for.
    
        if dim==1
          FinalTraj.Trajectories{i}=[FinalTraj.Trajectories{i}(:,1)*umperpx]; 
        elseif dim==2          
          FinalTraj.Trajectories{i}=[FinalTraj.Trajectories{i}(:,1)*umperpx FinalTraj.Trajectories{i}(:,2)*umperpx];
        elseif dim==3
          FinalTraj.Trajectories{i}=[FinalTraj.Trajectories{i}(:,1)*umperpx FinalTraj.Trajectories{i}(:,2)*umperpx FinalTraj.Trajectories{i}(:,3)*umperpx]; %ones(length(FinalTraj{i}(:,1)),1)];
        end

        if size(FinalTraj.Trajectories{i},1)<lengthlowerbound || size(FinalTraj.Trajectories{i},1)>lengthhigherbound
          FinalTraj.Trajectories{i}=[];
        end
        
end

%remove empty cells:

FinalTraj.Trajectories=FinalTraj.Trajectories(~cellfun('isempty',FinalTraj.Trajectories));

% Linear interpolation for NaN values in the matrices. The NaN values
% indicate the gaps of the fluorophore.
        
for i=1:size(FinalTraj.Trajectories,2)
        
        X=FinalTraj.Trajectories{i}(:,1);
        X=fixgaps(X);
        
        Y=FinalTraj.Trajectories{i}(:,2);
        Y=fixgaps(Y);
        
        FinalTraj.Trajectories{i}(:,1)=X;
        FinalTraj.Trajectories{i}(:,2)=Y;

end

F{k}=FinalTraj;

end

end

