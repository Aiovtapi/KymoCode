function [tracks] = LionToMSD(d,N_particles,N_time_steps)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

N_dim = 2; % 2D

tracks = cell(N_particles, 1);

% Time step between acquisition; fast acquisition!
dT = 0.05; % s,
time = (0 : N_time_steps-1)' * dT;

for i=1:N_particles
tracks{i}=[time d.x{i}(:,2) d.x{i}(:,4)];
end

end

