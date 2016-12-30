function [Traj] = LionStart(N)
%% Roy's implementation
% loading utrack trajectories
% utrack to VB3 format
% saving in new structure variable

% experiment numbers

%Lower upperbound for trajectories
LB=3;
UB=200;

%Traj: all trajectories from experiments N, from utrack to VB3 format
Traj=LionPrepare(N,LB,UB);

end

