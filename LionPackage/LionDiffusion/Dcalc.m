function [D] = Dcalc(Traj,dt)
%DCALC Calculate apparent diffusion from trajectories

for i=1:size(Traj,2)
    N=size(Traj{i},1);
    MSD(i)=0;
    for j=1:N-1
    MSD(i)=MSD(i)+(Traj{i}(j+1,1)-Traj{i}(j,1))^2+(Traj{i}(j+1,2)-Traj{i}(j,2))^2;
    end
    MSD(i)=MSD(i)/(N-1);
    D(i)=MSD(i)/(4*size(Traj{i},1)*dt);
end

end

