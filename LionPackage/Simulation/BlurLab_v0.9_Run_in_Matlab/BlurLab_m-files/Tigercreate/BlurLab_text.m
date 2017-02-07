%BlurLab Input Fluorescence Text File Maker
%Tristan Ursell, 2011
%
%BlurLab_text(X,Y,Z,I,F,L,fname)
%
%X Y and Z are the columns vectors specifying the positions of all
%fluorescent objects.  I specifies those objects' intensities.  F is the
%frame number in which sub-groups of objects will appear.  L is an optional
%label for each fluorescent object that must be the same from one frame to
%another and must increase sequentially.  'fname' is the name of the output
%text file.
%

function BlurLab_text(varargin)

if nargin==6
    X=varargin{1};
    Y=varargin{2};
    Z=varargin{3};
    I=varargin{4};
    F=varargin{5};
    fname=varargin{6};
    
    sx=size(X,2);
    sy=size(Y,2);
    sz=size(Z,2);
    si=size(I,2);
    sf=size(F,2);
    
    if or(or(or(or(sx~=1,sy~=1),sz~=1),si~=1),sf~=1)
        error('Inputs must be column vectors.')
    end
    
    Nx=length(X);
    Ny=length(Y);
    Nz=length(Z);
    Ni=length(I);
    Nf=length(F);
    
    if or(or(or(Ny~=Nx,Nz~=Nx),Nx~=Ni),Nx~=Nf)
        error('Input vectors must be the same length.')
    end
    
    fid=fopen(fname,'w+');
    disp(['Creating file: ' fname ' in ' pwd ', with ' num2str(Nx) 'lines.'])
    fprintf(fid,'%6.3f %6.3f %6.3f %6.3f %6f\n',[X';Y';Z';I';F']);
    fclose(fid);

elseif nargin==7
    X=varargin{1};
    Y=varargin{2};
    Z=varargin{3};
    I=varargin{4};
    F=varargin{5};
    L=varargin{6};
    
    fname=varargin{7};
    
    sx=size(X,2);
    sy=size(Y,2);
    sz=size(Z,2);
    si=size(I,2);
    sf=size(F,2);
    sl=size(L,2);
    
    if or(or(or(or(or(sx~=1,sy~=1),sz~=1),si~=1),sf~=1),sl~=1)
        error('Inputs must be column vectors.')
    end
    
    Nx=length(X);
    Ny=length(Y);
    Nz=length(Z);
    Ni=length(I);
    Nf=length(F);
    Nl=length(L);
    
    if or(or(or(or(Ny~=Nx,Nz~=Nx),Nx~=Ni),Nx~=Nf),Nx~=Nl)
        error('Input vectors must be the same length.')
    end
    
    fid=fopen(fname,'w+');
    disp(['Creating file: ' fname ' in ' pwd ', with ' num2str(Nx) ' lines.'])
    fprintf(fid,'%6.3f %6.3f %6.3f %6.3f %6f %6f\n',[X';Y';Z';I';F';L']);
    fclose(fid);
else
    error('Incorrect number of inputs.')
end