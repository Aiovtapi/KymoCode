function [] = Dplot(runinputfile,Title)
%DPLOT two dimensional diffusion plot from VB3 analysis input file.
% add VB3 to path!

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

for i=1:size(Traj,2)
    N=size(Traj{i},1);
    MSD(i)=0;
    for j=1:N-1
    MSD(i)=MSD(i)+(Traj{i}(j+1,1)-Traj{i}(j,1))^2+(Traj{i}(j+1,2)-Traj{i}(j,2))^2;
    end
    MSD(i)=MSD(i)/(N-1);
    D(i)=MSD(i)/(4*size(Traj{i},1)*dt);
end

%% plot
Nbins=100;
bins=linspace(0, 1, Nbins);
Dc=histc(D,bins);

fig1=figure(1);
bar(bins,Dc/(sum(Dc)))
set(gca,'fontsize',18)
axis([-0.005 0.5 0 0.1])
xlabel('Diffusion Coefficient (um2/s)')
ylabel('Probability (-)')
title(strcat({Title},'',{ 'Ntraj ='}, num2str(size(Wbest.T,2))));
saveas(gcf,'Figures/Dplot.fig')

end

