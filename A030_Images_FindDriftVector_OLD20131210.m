function dum=ChannelsImages_FindDriftVector
%This code loads a movie, allows the user to click a region with a blob and
%tracks a vector from it.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%NOTES: 
%'dipimage' should be installed and started ('dipstart') to make this
%code work
%Jacob Kerssemakers 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; 
%dipstart;
%expno='TEST';
%expno='CM_DnaN_37Deg_Series1002';
%expno='CM_DnaN-Dif-Gamma-ve-Position1_Series1';
expno='D005_CM_DnaN_Dif_37_deg_Rep_Check_002_C2'; %Dit is een D005 Experiment

initval=A001_Images_Set_Experiment(expno); %define your paths and files
initval.domovieload=1;  %default=1 (new analysis) %if not saved workspace before
initval.hROI=22;

%I) First, collect an image series. 
ImagesWorkspaceName=strcat(initval.basepath,initval.outname,'_Images',num2str(initval.maxfile),'.mat');
load(ImagesWorkspaceName,'aa');

imx=max(max(aa(:,:,1))); imn=min(min(aa(:,:,1)));
dipshow(aa(:,:,1),[imn imx]);   %show first im
[r,c]=ginput(1);     %pick a ROI
lor=r-initval.hROI;
hir=r+initval.hROI;
loc=c-initval.hROI;
hic=c+initval.hROI;
bb=double(dip_array((aa(lor:hir,loc:hic,:))));
kernel=squeeze(bb(:,:,1));
kernel=abs(kernel-mean(mean(kernel)));


imx=max(max(kernel)); imn=min(min(kernel));
dipshow(kernel,[imn imx]);   %show first im

[r,c,ls]=size(aa);

driftvector=zeros(ls,2);


%Here we first verify whether we need to create the folder
FolderExistence = exist(strcat(initval.basepath,initval.FiguresFolder,'DriftVector_ROI'));
if FolderExistence == 0
    mkdir(strcat(initval.basepath,initval.FiguresFolder,'DriftVector_ROI_FFT'))
end

%Here we first verify whether we need to create the folder
FolderExistence = exist(strcat(initval.basepath,initval.FiguresFolder,'DriftVector_CR'));
if FolderExistence == 0
    mkdir(strcat(initval.basepath,initval.FiguresFolder,'DriftVector_CR'))
end

fkernel=(fft2(kernel));
figure;
for i=1:ls
    ls-i
    im=bb(:,:,i);
    [x0,y0,x,y,prfx,prfy]=Track_Kernel(im,fkernel,i,strcat(initval.basepath,initval.FiguresFolder,'DriftVector_ROI_FFT/'),strcat(initval.basepath,initval.FiguresFolder,'DriftVector_CR/'));
    driftvector(i,:)=[x0 y0];
end
driftvector(:,1)=driftvector(:,1)-driftvector(1,1);
driftvector(:,2)=driftvector(:,2)-driftvector(1,2);

H=figure;
hold on;
plot(driftvector(:,1),'-bo','LineWidth',4,...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor','w',...
                'MarkerSize',6);
plot(driftvector(:,2),'-ro','LineWidth',4,...
                'MarkerEdgeColor','r',...
                'MarkerFaceColor','w',...
                'MarkerSize',6);
%set(gca, 'YTick', [0 - 0.02 0.03]);
set(gca, 'fontsize', 26, 'linewidth', 4, 'fontweight', 'bold');
set(gca,'TickLength',[0.02 0.02]);
ylabel('Position (px)', 'fontsize', 26, 'fontweight', 'bold');
xlabel('Frames (-)', 'fontsize', 26, 'fontweight', 'bold');
legend('X-drift', 'Y-drift');

%Here we first verify whether we need to create the folder
FolderExistence = exist(strcat(initval.basepath,initval.FiguresFolder));
if FolderExistence == 0
    mkdir(strcat(initval.basepath,initval.FiguresFolder))
end

%Writing the drift vector plot to the FIgures folder at the dataset
FigureName ='DriftVectorPlot';
initval.FigureToBeWritten=strcat(initval.basepath,initval.FiguresFolder,FigureName);
print(H, '-dpdf', '-r600',initval.FigureToBeWritten)
hold off;
drift=driftvector;
lbl=strcat(initval.basepath,initval.driftfile);

dlmwrite(lbl,driftvector);
initval.ImagesWorkspaceName=strcat(initval.basepath,initval.outname,'_Images',num2str(initval.maxfile),'.mat');
save(initval.ImagesWorkspaceName, 'drift', '/append');
disp('done');
end

function [x0,y0,x,y,prfx,prfy]=Track_Kernel(im,fkernel,Nr,DirNameROI,DirNameCR); 
 %cross-correlates image with template image
     im=abs(im-mean(mean(im))); [r,c]=size(im);
     
     cr=abs(fftshift(ifft2(fft2(im).*fkernel')));      %cross-correlation
    
     subplot(2,2,1); pcolor(im); colormap bone; shading flat; 
     subplot(2,2,2); pcolor(cr); colormap bone; shading flat; pause(0.02);
  
     
     FigureName_ROI ='ROI_Img';
     FigureName_CR = 'CR_Img';
     initval.FigureToBeWritten_ROI=strcat(DirNameROI,FigureName_ROI, num2str(Nr), '.tif');
     initval.FigureToBeWritten_CR=strcat(DirNameCR,FigureName_CR, num2str(Nr), '.tif');
     %fNameToWrite = ['orrected/RBall_Illum_Corrected' num2str(k) '.tif']; 
     
     AMIN = 0;
     AMAX = 65535;
    
     fr_ROI = double(im);
     imwrite(uint16(65535*mat2gray(fr_ROI,[AMIN AMAX])),initval.FigureToBeWritten_ROI,'Compression','none');

     %fr_CR = double(cr);
     fr_CR=cr;
     imwrite(uint16(65535*mat2gray(fr_CR)),initval.FigureToBeWritten_CR,'Compression','none');
     
     [val,x0]=max(max(cr)); [val,y0]=max(max(cr')); %maximum image-centered
    
     prfx=mean(cr)';         prfy=mean(cr')';          %averaged crosslines, handy for blobs
    % x=subpix_aroundzero(prfx)+x0; y=subpix_aroundzero(prfy)+y0;
    x=x0; y=y0;
end

function  x=subpix_aroundzero(prfx);
     xax=[-4:1:2]'; [val,cx]=max(prfx); c=length(cx);
     xa=mod(xax+cx,c)+1; prfx=prfx(xa);   %peak parabols with edge transfer
     prms=polyfit(xax,prfx,2); x=-prms(2)/(2*prms(1));
     
end




