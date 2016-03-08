function [nwx,nwy,Ispot,Ibackground_level,spotim_clipped,bckim,totim,Pixels]=LionMasker(spotim,x,y,ClipmaskR,GaussmaskW)
  %This function follows LL-G et al to get a masked spot.
                    %note that we do not resample with the clipping mask
                    %(which is bad!)
                    %ClipmaskR=5;
                    %GaussmaskW=4;
                    
                    [rr,cc]=size(spotim);
                    [XX,YY]=meshgrid(1:cc,1:rr);
                    II=0*XX+1;  %no weights
                    radpos=((XX-x).^2+(YY-y).^2).^0.5;
                    GaussMask=abs(1-exp(-radpos/(2*(GaussmaskW)).^2));
                    
                    sel=find(radpos<ClipmaskR); % evt. ClipmaskR op basis interpolatie
                    unsel=find(radpos>=ClipmaskR); 
                   
                    clippedGauss=GaussMask;
                    clippedGauss(unsel)=1;
                    
                    spotim_clipped=spotim; 
                    spotim_clipped(unsel)=0; 
                    outsideim=spotim;
                    outsideim(sel)=0;
                    
%                     if size(nonzeros(spotim_clipped(:)),1)<9
%                     Ispot=0;
%                     Ibackground_level=0;
%                     spotim_clipped=0;
%                     totim=[];
%                     bckim=0;
%                     nwx=0;
%                     nwy=0;
%                     else
                    Ibackground_level=mean(outsideim(unsel));
                    spotim_bc=spotim_clipped;
                    spotim_bc(sel)=spotim_bc(sel)-Ibackground_level;
                    spotim_bc(spotim_bc<0)=0;

                    Ispot=sum(spotim_bc(:));
                    Pixels=size(nonzeros(spotim_bc),1);
                    if isnan(Ispot)
                        Ispot=0;
                    end
                     
                    spotim_masked=spotim_clipped.*GaussMask;
                    bckim=spotim-spotim_bc;
                    totim=spotim_masked+bckim;
                        
                    clippedpoints=[XX(sel) YY(sel) spotim_masked(sel)];       
                    %[nwx,nwy,~,~]=JKD2_XY_calculate2Dmomentpoints(clippedpoints,1);
                    nwx=0; nwy=0;
                   
                    if 0
                    hold on
                    imagesc(spotim_bc); 
                    plot(nwx+0.5,nwy+0.5,'ro');
                    hold off
                    dum=1;
                    end
                    
%                     end
end

