% GaussCalcs2

space_units = '?m';
time_units = 's';

pixelsize=0.18;

%data source, from gausscalcs
OBJ=SdnN;

%Which experiment?
L=4;

N_particles = size(OBJ{L}.x,2);
N_time_steps = size(OBJ{L}.x{1},1);

% Time step between acquisition; fast acquisition!
dT = 0.01; % s,

MA=msdanalyzer(2,space_units,time_units);

tracks=LionToMSD(OBJ{L},N_particles,pixelsize,dT);

%add tracks to db
MA=MA.addAll(tracks);

%computeMSD
MA=MA.computeMSD;
t = (0 : N_time_steps-1)' * dT;
[T1, T2] = meshgrid(t, t);
all_delays = unique( abs(T1 - T2) );

fprintf('Found %d different delays.\n', numel(all_delays));
fprintf('For %d time-points, found %d different delays.\n', N_time_steps, size(MA.msd{1}, 1 ) );

%% plots
ax = gca; % current axes
ax.FontSize = 16;
set(ax,'LineWidth',2);
%     MA.plotMSD
 MA.plotTracks
% MA.labelPlotTracks;
% MA.plotMeanMSD(gca, true)


%% diffusion coefficient
[fo, gof] = MA.fitMeanMSD;

