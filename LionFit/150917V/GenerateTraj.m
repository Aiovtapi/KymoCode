%GaussFitSimEdit has to be run to load results.
%
%% Load results into structure: Data

clear
close
clc

%experimental pars:

exp='Roy_MM_Tus_dif';
lionval.channr=7;
lionval.viewchan='CFP';
lionval.viewbac=1;
lionval=LionDefine(exp,lionval);

timestep=0.02;                      %time per frame
D.um=1;                             %upper limit estimate diffusion (um^2/s)
D.px=D.um/(0.16^2);                 %in pixels

MSD=2*D.px*timestep;                %Mean Squared Displacement

AC=4*MSD;                           % Spot Appear Cost (Distance^2)
DAC=2*MSD;                          % Spot Disappear Cost 

% Cell range
Cell=1;

for N=Cell;
        thisbacfolder=strcat(lionval.bacstring{1},num2str(N,'%03.0f'));
        Data{Cell}=load(strcat(lionval.Mainfolder,'Results',lionval.OSslash,thisbacfolder));
end

%% 
Timelength=[];
for i=Cell
    Timelength=[Timelength size(Data{i}.x{1},1)];
end

Tmax=max(Timelength);
Ncells=length(Cell);

Ctotald=[];
Clinkd=cell(Ncells,Tmax);
Vectorx=cell(Ncells,Tmax);
Vectory=cell(Ncells,Tmax);
VectorX=cell(Ncells,Tmax);
VectorY=cell(Ncells,Tmax);
VectorI=cell(Ncells,Tmax);
Vectori=cell(Ncells,Tmax);
AuxiliaryMatrix=cell(Ncells,Tmax);


for i=Cell;
    
    I=size(Data{N}.x{1},1);
    Nspots=size(Data{N}.x,2);
    
    for t=1:I-1;
        
        for j=1:Nspots; % Rows at t
            for k=1:Nspots; % Columns at t+1
            
                if ~Data{i}.x{j}(t,2)==0 && ~Data{i}.x{k}(t+1,2)==0
                DistanceSqrd{i,t}(j,k)=(sqrt(Data{i}.x{j}(t,2).^2+Data{i}.x{j}(t,4).^2) - ...
                    sqrt(Data{i}.x{k}(t+1,2).^2+Data{i}.x{k}(t+1,4).^2)).^2;
                else
                DistanceSqrd{i,t}(j,k)=0;
                end
                
            end
                Vectorx{i,t}=[Vectorx{i,t} Data{i}.x{j}(t+1,2)]; %Vectors used for matrix product with the LAP solution
                Vectory{i,t}=[Vectory{i,t} Data{i}.x{j}(t+1,4)]; %This will transfer the position data to according spots
                Vectori{i,t}=[Vectori{i,t} Data{i}.x{j}(t+1,1)];
        end
        
                appearCost=ones(Nspots,1)*AC;                   %Spot Appear Cost
                disappearCost=ones(Nspots,1)*DAC;               %Spot Dissappear Cost

                aMatrix=diag(appearCost);                       %Diagmatrix
                dMatrix=diag(disappearCost);
                
                AuxiliaryMatrix{i,t}=DistanceSqrd{i,t}'; % given the lower cost of the matrix so that it doesn't influence final solution

                DistanceSqrd{i,t}(j+1:j+Nspots,1:k)=aMatrix; %lower left matrix
                DistanceSqrd{i,t}(1:j,k+1:k+Nspots)=dMatrix; %upper right matrix
                DistanceSqrd{i,t}(j+1:j+Nspots,k+1:k+Nspots)=AuxiliaryMatrix{i,t};
                
                Ctotald=[Ctotald;DistanceSqrd{i,t}(:)];
                
                DistanceSqrd{i,t}(DistanceSqrd{i,t}==0)=1000;
                
                [Amind{i,t},Costd{i,t}]=LionLAP(DistanceSqrd{i,t});

                %Reassign x-position correspondingly (and corresp
                %intensities)
                %First 'double' the size
                    
                Vectorx{i,t}=[Vectorx{i,t} Vectorx{i,t}];
                Vectory{i,t}=[Vectory{i,t} Vectory{i,t}];
                Vectori{i,t}=[Vectori{i,t} Vectori{i,t}];
                
                %Matrix Product
                VectorX{i,t}=Amind{i,t}*Vectorx{i,t}';
                VectorY{i,t}=Amind{i,t}*Vectory{i,t}';
                VectorI{i,t}=Amind{i,t}*Vectori{i,t}';
                
                %Correspond results to spots
                
                for j=1:Nspots
                                        
                    Data{i}.x{j}(t+1,2)=VectorX{i,t}(j);
                    Data{i}.x{j}(t+1,4)=VectorY{i,t}(j);
                    Data{i}.x{j}(t+1,1)=VectorI{i,t}(j);
                    
                end
        
    end
end




