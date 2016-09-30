function [] = PlotTrajs(Traj,D,Dthreshold)
%PLOTTRAJS Summary of this function goes here
%   Detailed explanation goes here
if nargin<2
    dt=0.02;
    Dthreshold=0.02;
    D=Dcalc(Traj,dt);
end

SizeT=size(Traj,2);
Tfast=cell(1,SizeT);
Tslow=cell(1,SizeT);

Ps=0.159; % PixelSize 0.159 micron per pixel

figure(1)
hold on
for i=1:SizeT   
    Traj{i}=Traj{i}./Ps;
    
    if D(i)>Dthreshold
        Tfast{i}=Traj{i};
    else
        Tslow{i}=Traj{i};
    end
    
end

for i=1:SizeT
    if ~isempty(Tfast{i});
    plot(Tfast{i}(:,1),Tfast{i}(:,2),'b');
    end
end

for i=1:SizeT   
    if ~isempty(Tslow{i});
    plot(Tslow{i}(:,1),Tslow{i}(:,2),'r');
    end
end

end

