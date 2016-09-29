    
function [X0,Y0,Z0,XYC,T0,I0,L0,D0,pts,npts,NewL] = ...
    Tig_Appearance(X0,Y0,Z0,XYC,T0,I0,L0,D0,npts,meanI,NewL,ini,TLx,Ly,Lz)
    
    % Get random number of spots to appear
    N_app = round(poissrnd(ini.Apois));

    % Get new values and save to matrices
    X0 = [X0; TLx*rand(N_app,1)];
    Y0 = [Y0; Ly*rand(N_app,1)];
    Z0 = [Z0; Lz*rand(N_app,1)];
    I0 = [I0; poissrnd(meanI,N_app,1)];
    T0 = [T0; ini.Framesec*rand(N_app,1)];
    
    L_new = (npts:npts+N_app-1)';
    D_new = logical(randi([0,1],N_app,1));
    XYC_new = zeros(N_app,2);
    
    % Get new directions
    for i_new = 1:N_app
        if D_new(i_new)
            Daci = 1;
        else
            Daci = 2;
        end
        XYC_new(i_new,:) = 2*(randi([0,1],1,2)-1/2) .* ...
            normrnd(ini.Dac(Daci),ini.SigmaD(Daci),1,2);
    end
    
    % Save values
    L0 = [L0; L_new];
    D0 = [D0; D_new];
    XYC = [XYC; XYC_new];
    
    npts = npts + N_app;

    pts = length(X0);
    NewL = [NewL; L_new];
end
