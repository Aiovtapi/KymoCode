function [mmask, nmask,bacpic,croppedimg] = Createbac(init,imageframe,thismesh,thisBbox,thisbacsize,bacpath,frami)
            
    bound = init.Extrabound;
    framesize = size(imageframe);
    thisRbox = thisBbox;
    Rbacsize = thisbacsize;
    thismsize = round(thisbacsize/init.pcresize);

    % Create mask from mesh
    thismeshl = [thismesh(:,1:2);thismesh(:,3:4)] + bound;
    thiscropmesh = round([thismeshl(:,1)-thisBbox(1), thismeshl(:,2)-thisBbox(3)]/init.pcresize);
    mask = poly2mask(thiscropmesh(:,1)',thiscropmesh(:,2)',thismsize(2),thismsize(1));
    
    % Dilate mask to get entire spot
    dimask = imdilate(mask,strel('disk',init.strelval));

    % Add bound to imageframe
    bound = floor(bound);
    nimgframe = padarray(imageframe,[bound,bound]);
    
    % Create non-masked bacpic from Cellbox
%     croppedimg = nimgframe(thisRbox(3):floor(thisRbox(3)+Rbacsize(2)-1),thisRbox(1):floor(thisRbox(1)+Rbacsize(1)-1));
    
    croppedimg = imcrop(nimgframe, [thisRbox(1), thisRbox(3), Rbacsize(1) - 1, Rbacsize(2) - 1]);

    % Ensure that size of the mask is the same as the size of the bacpic
    cimgsize = size(croppedimg);
    nmask = double(imresize(dimask,cimgsize));
    mmask = double(imresize(mask,cimgsize));
    
    % Creating the masked bacpic by applying the mask to the bacpic
    bacpic = uint16(nmask.*croppedimg);
    thisbacpath = strcat(num2str(frami,'%03.0f'),'.tif');
    imwrite(bacpic,strcat(bacpath,thisbacpath));
    


%     % Add bound to imageframe
%     bound = floor(bound);
%     nimgframe = padarray(imageframe,[bound,bound]);
%     
%     % Create non-masked bacpic from Cellbox
%     croppedimg = nimgframe(thisRbox(3):floor(thisRbox(3)+Rbacsize(2)-1),thisRbox(1):floor(thisRbox(1)+Rbacsize(1)-1));
% 
%     % Ensure that size of the mask is the same as the size of the bacpic
%     cimgsize = size(croppedimg);
%     
%     % Create mask from mesh
%     thismeshl = [thismesh(:,1:2);thismesh(:,3:4)] + bound;
%     thiscropmesh = round([thismeshl(:,1)-thisBbox(1), thismeshl(:,2)-thisBbox(3)]/init.pcresize);
%     mmask = poly2mask(thiscropmesh(:,1)',thiscropmesh(:,2)',cimgsize(1),cimgsize(2));
%     
%     % Dilate mask to get entire spot
%     nmask = imdilate(mmask,strel('disk',init.strelval));
    

end