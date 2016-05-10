clear all

dif1=load('1/RBIC/results/difresults.mat');
dif2=load('2/RBIC/results/difresults.mat');

Tus1=load('1/RBIC/results/tusresults.mat');
Tus2=load('2/RBIC/results/tusresults.mat');

DnaN1=load('2/RBIC/results/dnanresults.mat');

Ncells1=size(dif1.cellList.meshData{1},2);
Ncells2=size(dif2.cellList.meshData{1},2);

notusspots=0;
I=0;
for i=1:Ncells1
    I=I+1;
    if isfield(dif1.cellList.meshData{1}{i},'spots') && isfield(Tus1.cellList.meshData{1}{i},'spots')
    SpotsVardif{i}=dif1.cellList.meshData{1}{i}.spots.l;
    SpotsVarTus{i}=Tus1.cellList.meshData{1}{i}.spots.l;
    TusSpotNumber(i)=size(Tus1.cellList.meshData{1}{i}.spots.x,2);
    difSpotNumber(i)=size(dif1.cellList.meshData{1}{i}.spots.x,2);
    MeshLength(i)=length(dif1.cellList.meshData{1}{i}.mesh(:,1));
    [CellLength(i),~]=projectToMesh(dif1.cellList.meshData{1}{i}.mesh(MeshLength(i),1),dif1.cellList.meshData{1}{i}.mesh(MeshLength(i),2),...
        dif1.cellList.meshData{1}{i}.mesh);
    NormTusPos{i}=SpotsVarTus{i}/CellLength(i);
    NormdifPos{i}=SpotsVardif{i}/CellLength(i);
    Indx(I)=i;
    else
        notusspots=notusspots+~isfield(Tus1.cellList.meshData{1}{i},'spots');
        continue
    end
end
for i=Ncells1+1:Ncells2
    I=I+1;
    if isfield(dif2.cellList.meshData{1}{i-Ncells1},'spots') && isfield(Tus2.cellList.meshData{1}{i-Ncells1},'spots')
    SpotsVardif{i}=dif2.cellList.meshData{1}{i-Ncells1}.spots.l;
    SpotsVarTus{i}=Tus2.cellList.meshData{1}{i-Ncells1}.spots.l;
    TusSpotNumber(i)=size(Tus2.cellList.meshData{1}{i-Ncells1}.spots.x,2);
    difSpotNumber(i)=size(dif2.cellList.meshData{1}{i-Ncells1}.spots.x,2);
    MeshLength(i)=length(dif2.cellList.meshData{1}{i-Ncells1}.mesh(:,1));
    [CellLength(i),~]=projectToMesh(dif2.cellList.meshData{1}{i-Ncells1}.mesh(MeshLength(i),1),dif2.cellList.meshData{1}{i-Ncells1}.mesh(MeshLength(i),2),...
    dif2.cellList.meshData{1}{i-Ncells1}.mesh);
    NormTusPos{i}=SpotsVarTus{i}/CellLength(i);
    NormdifPos{i}=SpotsVardif{i}/CellLength(i);
    Indx(I)=i;
    else
        notusspots=notusspots+~isfield(Tus2.cellList.meshData{1}{i-Ncells1},'spots');
        continue
    end
end

Indx=nonzeros(Indx)';

J=0;
NormTusPosMat=zeros(length(Indx),6);

quadtus=zeros(length(Indx),4);

for i=Indx
    if TusSpotNumber(i)==1;
        singletus(i,1)=NormTusPos{i};
    elseif TusSpotNumber(i)==2;
        doubletus(i,:)=NormTusPos{i}; singletus(i,:)=0;
    elseif TusSpotNumber(i)==3;
        tripletus(i,:)=NormTusPos{i}; singletus(i,:)=0; doubletus(i,:)=[0 0];
    elseif TusSpotNumber(i)==4;
        quadtus(i,:)=NormTusPos{i}; tripletus(i,:)=0; singletus(i,:)=0; doubletus(i,:)=0;
    elseif TusSpotNumber(i)==5;
        quinttus(i,:)=NormTusPos{i}; quadtus(i,:)=0; tripletus(i,:)=0; doubletus(i,:)=0; singletus(i,:)=0;
    elseif isempty(TusSpotNumber(i))
        continue
    end
end

for i=Indx
    if difSpotNumber(i)==1;
        singledif(i,1)=NormdifPos{i};
    elseif difSpotNumber(i)==2;
        doubledif(i,:)=NormdifPos{i};
    elseif difSpotNumber(i)==3;
        tripledif(i,:)=NormdifPos{i};
    elseif isempty(difSpotNumber(i))
        continue
    end
end

% l : coordinate centerline
% d : distance from centerline
% x : euclidean coordinate in image
% y : euclidean coordinate in image
% position: segment number in which spot is located (if 0 outside of cell)
% adj.Rsquared: goodness of fit, 1 indicates almost perfect.
% confidenceInterval_x_y : 95% confidence bounds of parameters c and y of
% fit. x and y being coordinate values of spot.

