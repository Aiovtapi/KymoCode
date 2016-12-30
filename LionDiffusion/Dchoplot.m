function [D,TrajC] = Dchoplot(runinputfile,Title)
%DCHOPLOT Summary of this function goes here
%   Detailed explanation goes here

if nargin<2
    Title='';
end

if(ischar(runinputfile))
    if(exist(runinputfile,'file')==2)
        opt=VB3_getOptions(runinputfile);
        %disp(['Read runinput file ' runinputfile])
    else 
       error(['Cannot find runinput file ' runinputfile '. Please check name and path.'])
    end
elseif(isstruct(runinputfile))
    opt=runinputfile;
    runinputfile=opt.runinputfile;
    %disp(['Read options structure based on runinput file ' runinputfile ])
else
    error('Could not find runinputfile or interpret runinputfile argument.')
end
filename=opt.outputfile;
load(filename);

Traj=importdata(opt.inputfile);
dt=opt.timestep;
TrajC=cell(1,size(Traj,2));

for i=1:size(Traj,2)
    N=size(Traj{i},1); %trajectory length
    if N>5
    Rem=mod(N,5); %remainder of traj after chopping
    Nloop=(N-Rem)/5; 
    NC=size(TrajC,2);
    TrajC{i}=Traj{i}(1:5,1:2);
    for j=2:Nloop
        TrajC{j-1+NC}=Traj{i}(j*5-4:j*5,1:2);
    end
    TrajC{NC+Nloop}=Traj{i}(Nloop*5+1:Nloop*5+Rem,1:2);
    else
        TrajC{i}=Traj{i};
    end
end

TrajC=TrajC(~cellfun('isempty',TrajC));

for i=1:size(TrajC,2);
    Nc=size(TrajC{i},1);
    MSD(i)=0;
    if Nc>3
    for j=1:Nc-1
    MSD(i)=MSD(i)+(TrajC{i}(j+1,1)-TrajC{i}(j,1))^2+(TrajC{i}(j+1,2)-TrajC{i}(j,2))^2;
    end
    MSD(i)=MSD(i)/(Nc-1);
    D(i)=MSD(i)/(4*Nc*dt);
    else
        TrajC{i}=[];
    end
end
MSD=nonzeros(MSD);
D=nonzeros(D);
D=D(D>0.015);

%% plot
Nbins=100;
bins=linspace(0, 2, Nbins);
Dc=histc(D,bins);

fig1=figure(1);
bar(bins,Dc/(sum(Dc)))
set(gca,'fontsize',18)
axis([-0.005 1.5 0 0.12])
xlabel('Diffusion Coefficient (um2/s)')
ylabel('Probability (-)')
title(strcat({Title},'',{ 'Ntraj ='}, num2str(size(D,1))));
saveas(gcf,'Figures/Dplot.fig')

end

