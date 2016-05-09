SPACE_UNITS = '?m';
TIME_UNITS = 's';

N_PARTICLES = 10;
N_TIME_STEPS = 100;
N_DIM = 2; % 2D

% Typical values taken from studies of proteins diffusing in membranes:
% Diffusion coefficient
D  = 1e-3; % ?m^2/s
% Time step between acquisition; fast acquisition!
dT = 0.05; % s,

% Area size, just used to disperse particles in 2D. Has no impact on
% analysis.
SIZE = 2; % ?m

k = sqrt(2 * D * dT);

tracks = cell(N_PARTICLES, 1);

for i = 1 : N_PARTICLES

    % Time
    time = (0 : N_TIME_STEPS-1)' * dT;

    % Initial position
    X0 = SIZE .* rand(1, N_DIM);

    % Integrate uncorrelated displacement
    dX = k * randn(N_TIME_STEPS, N_DIM);
    dX(1, :) = X0;
    X = cumsum(dX, 1);

    % Store
    tracks{i} = [time X];

end
clear i X dX time X0

ma=msdanalyzer(2,SPACE_UNITS,TIME_UNITS);

ma = ma.addAll(tracks);

% Compute MSD
ma = ma.computeMSD;
ma.msd

% Delays
t = (0 : N_TIME_STEPS)' * dT;
[T1, T2] = meshgrid(t, t);
all_delays = unique( abs(T1 - T2) );

fprintf('Found %d different delays.\n', numel(all_delays));
fprintf('For %d time-points, found %d different delays.\n', N_TIME_STEPS, size( ma.msd{1}, 1 ) );

% diffusion coefficient
[fo, gof] = ma.fitMeanMSD;

