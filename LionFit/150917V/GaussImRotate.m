function GaussImRotate(bacseriepth)

D=readtimeseries(strcat(bacseriepth,'.tif'),'tif');
data=dip_array(D);

for i=1:size(data,3)
    rdata(:,:,i)=imrotate(data(:,:,i),180);
    imwrite(rdata(:,:,i),strcat(bacseriepth,num2str(i,'%03.0f'),'.tif'),'tif');
end

disp('Bac series flipped')

end