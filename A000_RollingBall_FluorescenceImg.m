function A000_RollingBall_FluorescenceImg(expName)
%%M.C. Moolman, 2013
%%-This function performs a type of Rolling Ball filtering. It performs an
%%opening of the image, namely an erosion followed by a dilation using a
%%ball of defined radius as the structural element

%%Literature reference for a Rolling ball filter: 
%%S. Sternberg,Biomedical Image Processing, IEEE Computer, January 1983. 

%clc;
%clear all;
initval = A001_Set_Experiment_PAmCherry(expName);

%%Here we load the ball radius/height and the structural element
BallHeight = initval.BallRadius;
se = strel('ball', initval.BallRadius, BallHeight, 0);


fNameFL_Org= strcat(initval.BasePath,initval.FL_File_Original);
fNameRollingBalledtoWrite = strcat(initval.BasePath,initval.FL_Path_RBolled,initval.FL_File_RBolled_File)
AMIN = 0;
AMAX = 65535;

num_images = (initval.maxfile_Original - initval.minfile_Original);


FolderExistence = exist(strcat(initval.BasePath,initval.FL_Path_RBolled));
if FolderExistence == 0
    mkdir(strcat(initval.BasePath,initval.FL_Path_RBolled));
end

%%Here we read in the individual images, performs the Rolling ball
%%filtering and then write it to a tiff stack
for k = initval.minfile_Original:initval.maxfile_Original
strcat('Processing: frame ',' ',int2str(k-initval.minfile_Original+1),' of',int2str(num_images+1))
%for k = 13:18  
    SubStringNr = sprintf('%03i',k);
    FileToLoad = strrep(fNameFL_Org, '*', SubStringNr);
    I = imread(FileToLoad);
    
%     figure;
%     imagesc(I)
%     axis equal tight
%     box on
%     colormap(gray)
%     axis xy;
    
    BackgndImg = imopen(I,se);
    
    RollingBalled_Img = imsubtract(I,BackgndImg);
%     figure;
%     imagesc(RollingBalled_Img)
%     axis equal tight
%     box on
%     colormap(gray)
%     axis xy;

    %ImagesToWrite(:,:,k) = RollingBalled_Img;
    imwrite(uint16(65535*mat2gray(RollingBalled_Img,[AMIN AMAX])),fNameRollingBalledtoWrite,'WriteMode','append','Compression','none');
    
end