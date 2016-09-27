    
function [X0,Y0,Z0,I0,L0,A0,R0,D0,pts,npts,NewL] = ...
    Tig_Appearance(X0,Y0,Z0,I0,L0,A0,R0,D0,npts,meanI,NewL,Apois,Lx,Ly,Lz,D)

    N_app = round(poissrnd(Apois));

    L_new = (npts:npts+N_app-1)';

    X0 = [X0; Lx*rand(N_app,1)];
    Y0 = [Y0; Ly*rand(N_app,1)];
    Z0 = [Z0; Lz*rand(N_app,1)];
    I0 = [I0; poissrnd(meanI,N_app,1)];
    L0 = [L0; L_new];
    A0 = [A0; 2*pi*rand(N_app,1)];
    if N_app > 0
        R0 = [R0; normrnd(0,sqrt(2*D),N_app,1)];
    end
    D0 = [D0, logical(randi([0,1],1,N_app))];
    
    npts = npts + N_app;

    pts = length(X0);
    NewL = [NewL; L_new];
end
