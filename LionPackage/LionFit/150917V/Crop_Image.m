function [im,borders] =Crop_Image(im0)
        %Crop to nonzero area
            [r,c]=size(im0);
            [X,Y]=meshgrid(1:c,1:r);
            sel=find(im0~=0);
            lor=min(Y(sel));
            hir=max(Y(sel));
            loc=min(X(sel));
            hic=max(X(sel));
            im=im0(lor:hir,loc:hic);
            borders=[lor,hir,loc,hic];
