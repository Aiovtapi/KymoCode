function RoicotrascaPH(Imgspath,Imgname,Resizeval,Transval)
%% Presets
if nargin < 2;
    Resizeval = 1;
    Transval = [0,0];
end

%%

flpth = strcat(Imgspath,Imgname);
FLinfo = imfinfo(flpth);
num_images = numel(FLinfo);
Outpath = strcat(Imgspath,'Transformed_',Imgname);


for k = 1:num_images;
    disp(['Processing image ',num2str(k),' out of ',num2str(num_images)])
    
    FLimg = im2double(imread(flpth,k,'Info',FLinfo));
    
    %% Transformation
    
    Timg = imtranslate(imresize(FLimg,Resizeval),Transval,'FIllvalues',0);
    
    %% Write to file
    imwrite(uint16(65535*Timg),Outpath,'WriteMode','append','Compression','none');

end

disp('Done')
end
