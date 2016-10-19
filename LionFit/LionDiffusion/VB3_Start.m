%% Roy's implementation
% loading utrack trajectories
% utrack to VB3 format
% saving in new structure variable

% experiment numbers
N=[3 5];

%Lower upperbound for trajectories
LB=4;
UB=60;

%Traj: all trajectories from experiments N, from utrack to VB3 format
Traj=LionPrepare(N,LB,UB);