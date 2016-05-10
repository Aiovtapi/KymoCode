%% 2D Gaussian fit to spot -- Gaussplosion!!
%
% This code is written to identify multiple spots within the cell. The
% input are the Bacpic images, as a result of the KymoCode, of seperate 
% cells over time. 
%
% written by Roy de Leeuw - 2015
%
% output : x{j}(i,:) = [A,x,sx,y,sy,0]
%
% j : The spot number. Each image could have multiple spots. Stored in cell
% i : The image number. Each cell is image stack as function of time.
% A : Amplitude of the spot.
% x : x-coordinate of the spot. with x = 0 being the cell pole.
% sx : std in x of the spot.
% y : y-coordinate of the spot. 
% sy : std in y of the spot.
%
%% 

clear all
close all
clc
tic

%% Inputs 

exp='Roy_MM_Tus_dif';

lionval.channr=7;
lionval.viewchan='CFP';
lionval.viewbac=1:10;

lionval=LionDefine(exp,lionval);

    
for Cell=lionval.viewbac;
    
    imgflip = 1;        % If you want to align the CFP images.
    imgflipped = 0;
    
    while imgflip == 1; 

        clear XNorm x NSpots SNR ydatacrpd pixels

        ClipFactor=1;
        GaussFactor=1;

        thisbacfolder=strcat(lionval.bacstring{1},lionval.bacstring{2},num2str(Cell,'%03.0f'));
        bacseriepth=strcat(lionval.Mainfolder,thisbacfolder,lionval.OSslash);
        
        % Check for whether flipping is needed based on the difchannel for
        % non difchannels. 
        if ~strcmp(lionval.viewchan,lionval.difchan)
            if exist(strcat(lionval.diffolder,'Results',lionval.OSslash,thisbacfolder,'.mat'))
                flipimage=load(strcat(lionval.diffolder,'Results',lionval.OSslash,thisbacfolder,'.mat'),'imgflipped');
                if flipimage.imgflipped==1
                    disp(['Flipping ',thisbacfolder])
                    GaussImRotate(bacseriepth);  
                end
            else
                disp('Dif-channel not computed. Cell orientation not defined. Do you want to continue? Click to continue')
                waitforbuttonpress; close all
            end
            imgflip = 0;
            imgflipped = flipimage.imgflipped;           
        end
        
        
        d1{1}=readtimeseries(strcat(bacseriepth,'.tif'),'tif');

        data=dip_array(d1{1}); %turn into uint16 array

        %% Ze Defs

        Nspots=1; %ignore this name

        Tsize=size(data,3);

        XSize=zeros(size(data,2),1);
        YSize=zeros(size(data,1),1); %Ysize
        ydata=cell(Tsize,1); 
        Ydata=cell(Tsize,1); Ydata_X=zeros(Nspots,1); Ydata_Y=zeros(Nspots,1);
        ydatacrpd=cell(Tsize,1); ydatacrpdR1=cell(Tsize,1); 
        Ampguess=zeros(Tsize,1); Case=cell(Nspots,1); 
        x0=cell(Nspots,1); x=cell(Nspots,1); X=cell(Nspots,1); 
        Y=cell(Nspots,1); Yg=zeros(Tsize,Nspots); Xg=zeros(Tsize,Nspots); 
        xdata=cell(2,Nspots);
        Size=cell(Tsize,1); pixels=cell(5,1);

        lb = [0,0,0,0,0]; %lower bounds for fit

        % fit using levenberg-marquardt algorithm
        OPTIONS = optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','Display','off');

        %% Load Data
        % parameters : [Amplitude, x0, sigmax, y0, sigmay, angel(in rad)]

        % Define some parameters 

        SA=3; % pixel step from maximum amplitude 
        Sx=3; Sy=3; % guess of sigmaX of spots
        Px=3; Py=3; % for optimal fitting add black 'pads' around detected spot
        Bs=1; Bss=2*Bs+1; % this defines size of black sqaure
        lob=1; upb=size(data(:,:,1),1); chanthickness=size(data(:,:,1),1); % define the boundaries of the channel


        % Initiate the spotting process..

        i=1;

        %Determine outliers for determining intensity threshold.

        lowerboundchannel=7;


        Tolerance=2;
        sigmachange=0.5;
        how='positive';
        show=0;


        for i=1:Tsize

            %Remove the x,y-Pads
            ydata{i}=double(data(:,:,i));

        %   Noticed that cropping causes shift in simulations!!
        if lionval.cropindx==1;
           [ydatacrpd{i},~]=Crop_Image(ydata{i});
        else 
            ydatacrpd{i}=ydata{i};
        end
            ydatacrpdR1{i,1}=ydatacrpd{i};

            higherboundchannel=16;
            
            % Intensity thresholding for outliers
            Data=ydatacrpd{i}(lowerboundchannel:higherboundchannel,:);
            [flag,A]=DetermineOutliers(Data,Tolerance,sigmachange,how,show);
            Outliersdata=~flag.*ydatacrpd{i}(lowerboundchannel:higherboundchannel,:);
            IntensityPeakThreshold = mean(nonzeros(Outliersdata(:)))+std(nonzeros(Outliersdata(:)));

            ydatacrpdR1{i,1}=Outliersdata;
            
            j=1;


            [x0{j}(i,:),Case{j}(i),ydatacrpdR1{i,j+1},Ydata{i,j},Size{i,j},Yg(i,j),Xg(i,j)]= ... 
                LionSpotter(ydatacrpdR1{i,j},SA,Sx,Sy,Px,Py,Bs,lob,upb);
            
           %still needs background level to be mediated from within the channel.   
            while x0{j}(i,1)>IntensityPeakThreshold && j<20 && Size{i,j}(1)>0

            % make guesses for spot position

            % output = [guessvector,case,newcutoutimage,imageusedforfitting]
            % guessvector = [ampl,x,sigx,y,sigy]
            % case = 2 (spot out of channel), -1 (left edge), 0 (middle), 1 (right
            % edge)
            % imageusedforfitting = only the spot is used for fitting (cut out)

            % preparing framework for gaussian fitting
            [Ydata_X(j),Ydata_Y(j)]=size(Ydata{i,j});
            [X{j},Y{j}] =  meshgrid(linspace(1,Ydata_Y(j),Ydata_Y(j)), ...
                linspace(1,Ydata_X(j),Ydata_X(j)));
            xdata{1,j} = X{j};
            xdata{2,j} = Y{j};

            % define upper boundary (Size{i,j}(2) is the X size of the crpd-image)
            ub = [realmax('double'),Size{i,j}(2),(Size{i,j}(2))^2,Size{i,j}(2), ...
                (Size{i,j}(2))^2]; %upper bounds

            % Do the Gaussian fitting

            [x{j}(i,:),resnorm,residual,exitflag] = lsqcurvefit(@GaussPlosFunc, ...
               x0{j}(i,:),xdata(:,j),Ydata{i,j},lb,ub,OPTIONS);


           % X-SA-1 is the length in X to the point X-SA on the image.

            if Case{j}(i)==2
                x{j}(i,2)=x{j}(i,2)-Px;
                x{j}(i,4)=x{j}(i,4)-Py+(Yg(i,j)-SA)-1;
            elseif Case{j}(i)==3
                x{j}(i,2)=x{j}(i,2)-Px;
                x{j}(i,4)=x{j}(i,4)-Py+(Size{i,j}(1)-2*SA)-1;
            elseif Case{j}(i)==1
                x{j}(i,2)=x{j}(i,2)-Px;
                x{j}(i,4)=x{j}(i,4)-Py;
            elseif Case{j}(i)==9
                x{j}(i,2)=x{j}(i,2)-Px+(Xg(i,j)-SA)-1;
                x{j}(i,4)=x{j}(i,4)-Py+(Yg(i,j)-SA)-1;
            elseif Case{j}(i)==6
                x{j}(i,2)=x{j}(i,2)-Px+(Size{i,j}(2)-2*SA)-1;
                x{j}(i,4)=x{j}(i,4)-Py+(Yg(i,j)-SA)-1;
            elseif Case{j}(i)==7
                x{j}(i,2)=x{j}(i,2)-Px+(Size{i,j}(2)-2*SA)-1;
                x{j}(i,4)=x{j}(i,4)-Py;
            elseif Case{j}(i)==8
                x{j}(i,2)=x{j}(i,2)-Px+(Xg(i,j)-SA)-1; 
                x{j}(i,4)=x{j}(i,4)-Py;
            elseif Case{j}(i)==4
                x{j}(i,2)=x{j}(i,2)-Px+(Xg(i,j)-SA)-1;
                x{j}(i,4)=x{j}(i,4)-Py-1+(Size{i,j}(1)-2*SA);
            elseif Case{j}(i)==5
                x{j}(i,2)=x{j}(i,2)-Px+(Size{i,j}(2)-2*SA)-1;
                x{j}(i,4)=x{j}(i,4)-Py-1+(Size{i,j}(1)-2*SA);        
            elseif Case{j}(i)==10
                x{j}(i,1:8)=0;
                XNorm{j}(i,1:8)=0;
            end
            
           % Check for whether flipping is necessary. This will be skipped
           % for non dif-channels, as flipping has already been done. 
           if imgflip==1
                imgflip = x{1}(1,2)>size(ydatacrpd{1},2)/2;
           end
           
           % Prevent flipping loops. Only one round of flipping is allowed.
           if imgflipped==1
               imgflip = 0;
           end
           
           if imgflip == 1
               imgflipped = 1;
               break
           end



        %   Size{i,j}(1) is the height of the image.
        %   Size{i,j}(2) is the width of the image.

            XNorm{j}(i,2)=x{j}(i,2)/Size{i,j}(2);

            XNorm{j}(i,4)=(x{j}(i,4))/Size{i,j}(1);

            GaussmaskW=GaussFactor*(x{j}(i,3).^2+x{j}(i,5).^2)^(1/2);
            ClipmaskR=ClipFactor*(x{j}(i,3).^2+x{j}(i,5).^2)^(1/2);

           [~,~,ISPOT,Ibacklevel,spotim_clipped,bckim,ydatacrpdR1{i,j+1},pixels{j}(i)]=LionMasker(ydatacrpdR1{i,j},x{j}(i,2),x{j}(i,4),ClipmaskR,GaussmaskW);

           if size(nonzeros(spotim_clipped(:)),1)<1
               break
           else
            j=j+1;

                [x0{j}(i,:),Case{j}(i),~,Ydata{i,j},Size{i,j},Yg(i,j),Xg(i,j)]= ... 
                LionSpotter(ydatacrpdR1{i,j},SA,Sx,Sy,Px,Py,Bs,lob,upb);  
           end

            end

            % Breaking the loop for image flipping. 
            if imgflip==1
                break
            end
        end
        
        % Flipping the images
        if imgflip==1
            disp(['Flipping ',thisbacfolder])
            GaussImRotate(bacseriepth);
        end
    end

    NSpots=size(x,2);
    LXNormOne=size(XNorm{1},1);
    LxOne=size(x{1},1);

    % make sure that each spot data is same length
    % because some spots will be there for just a few frames
    % and their last will be the length of the array.

    for j=2:NSpots;
        x{j}=[x{j};zeros(LxOne-size(x{j},1),5)];
        XNorm{j}=[XNorm{j};zeros(LXNormOne-size(XNorm{j},1),4)];
    end


    %% Translate back to original coords

    % use XNorm for normalization to position w.r.t. cell

    %% Integrated Intensity V2.0
    II=cell(5,5,Tsize);
    FCII=cell(1,Tsize);
    ROI=cell(Tsize,Nspots);
    Ispot=zeros(Tsize,Nspots);

    %SNR vars
    SNR=zeros(Tsize,Nspots);
    SNRSimonetti=zeros(Tsize,Nspots);
    SNRBlur=zeros(Tsize,Nspots); %theoretical SNR of noise simulator

    %BlurNoise Parameters
    MeanNoise=1500;
    STDNoise=225;

    Gain=1;


    for i=1:Tsize

         % Full Cell Integrated Intensity
         FCII{i}=ydatacrpd{i}(lowerboundchannel:higherboundchannel,:);

        for j=1:NSpots

            %This is the full cell integrated intensity
            x{j}(i,7)=sum(sum(FCII{i}));
            % to do: compare spot positions
            if ~isempty(Size{i,j}) % there still has to be an image.

            padx=5;
            pady=5;

            [r,c]=size(ydatacrpdR1{i,j});
            [xx_ori,yy_ori]=meshgrid(-5:5,-5:5);


            GaussmaskW=GaussFactor*(x{j}(i,3).^2+x{j}(i,5).^2)^(1/2);
            ClipmaskR=ClipFactor*(x{j}(i,3).^2+x{j}(i,5).^2)^(1/2);

            ydatacrpdR1{i,j}=[zeros(r,padx) ydatacrpdR1{i,j} zeros(r,padx)]; %make pads because spots can be on the edge of the image
            ydatacrpdR1{i,j}=[zeros(pady,size(ydatacrpdR1{i,j},2)); ...
            ydatacrpdR1{i,j};zeros(pady,size(ydatacrpdR1{i,j},2))];
            [R,C]=size(ydatacrpdR1{i,j});
            [XX,YY]=meshgrid(1:C,1:R);

    %         if x{j}(i,2)>3 && x{j}(i,2)<Size{i,j}(2)-3
    %         xx0=x{j}(i,2)+xx_ori;
    %         elseif x{j}(i,2)>Size{i,j}(2)-3 && Size{i,j}(2)>1
    %         xx0=Size{i,j}(2)-3+xx_ori;
    %         else
    %         xx0=xx_ori+3;
    %         end

            xx0=x{j}(i,2)+padx+xx_ori;
            yy0=x{j}(i,4)+pady+yy_ori;

            ROI{i,j}=interp2(XX,YY,ydatacrpdR1{i,j},xx0,yy0);

            [rROI,cROI]=size(ROI{i,j});

            xROI=cROI/2;
            yROI=rROI/2;

            [~,~,Ispot(i,j),Ibackground_level(i,j),spotim_masked,bckim]=DoubleMaskedCom(ROI{i,j},xROI...
                ,yROI,ClipmaskR,GaussmaskW);

            RadiusSpot(i,j)=GaussmaskW;

            nn=size(spotim_masked(spotim_masked>0));

            BCKIM=bckim(1:5,:);

            SNR(i,j)=mean(spotim_masked(spotim_masked>0))/std(BCKIM(:));

            %SNRSimonetti(i,j)=(sum(spotim_masked(spotim_masked>0))-nn(1)*Ibackground_level)...
            %    /sqrt((mean(spotim_masked(spotim_masked>0))-nn(1)*Ibackground_level)/Gain)+nn(1)*std(bckim(:))^2+(nn(1)*std(bckim(:))^2)/size(bckim(:),1);

            %SNRBlur(i,j)=(sum(spotim_masked(spotim_masked>0))-nn(1)*MeanNoise)...
            %    /sqrt((mean(spotim_masked(spotim_masked>0))--nn(1)*MeanNoise)/Gain)+nn(1)*STDNoise^2+(nn(1)*STDNoise^2)/size(bckim(:),1);

            else 

            ROI{i,j}=0; SNR(i,j)=0;
            nwx{i,j}=0; nwy{i,j}=0; Ispot(i,j)=0; Ibackground_level(i,j)=0;
            spotim_masked=0;
            bckim=0;

            end

            x{j}(i,8)=Ispot(i,j);

        end
    end

    %% Integrated Intensity Calculation


    % for i=1:Zsize
    %     
    %     % Full Cell Integrated Intensity
    %     FCII{i}=ydatacrpd{i}(lob:upb,1:Size{i,j}(2));
    %     
    %     for j=1:Nspots
    %     
    %     %This is the full cell integrated intensity
    %     x{j}(i,7)=sum(sum(FCII{i}));
    %     
    %         % define bordercases
    %         % when x it is too close to the left border of image
    %         if round(x{j}(i,2))<=1
    %         BorderCase=-1;
    %         % when x is too close to the right border of image
    %         elseif round(x{j}(i,2))>=Size{i,j}(2)-1 
    %         BorderCase=1;
    %         else 
    %         BorderCase=0;
    %         end
    %     
    %         
    %         %define bordercases for wider spot
    %         if round(x{j}(i,2))<=3
    %         BorderCaseWide=-1;
    %         % when x is too close to the right border of image
    %         elseif round(x{j}(i,2))>=Size{i,j}(2)-3
    %         BorderCaseWide=1;
    %         else 
    %         BorderCaseWide=0;
    %         end
    %     % SPOT around centroid
    %     % If sigmaX == 2 or sigmaX == 1
    %     
    %     %For BorderCase == 0
    %     if (round(x{j}(i,3))==2 || round(x{j}(i,3))==1) && BorderCase ==0;
    %     
    %         for n=[-1 0 1] 
    %             if round(x{j}(i,5))==2 || round(x{j}(i,5))==1
    %                 for k=[-1 0 1]
    %                     II{n+2,k+2,i}=ydatacrpd{i}(round(x{j}(i,4))+k,round(x{j}(i,2))+n);
    %                 end
    %             end
    %             if round(x{j}(i,5))==3 || round(x{j}(i,5))==4
    %                 for k=[-2 -1 0 1 2]
    %                     II{n+2,k+3,i}=ydatacrpd{i}(round(x{j}(i,4))+k,round(x{j}(i,2))+n);
    %                 end
    %             end
    %         end
    %     
    %     % For BorderCase == -1
    %     elseif (round(x{j}(i,3))==2 || round(x{j}(i,3))==1) && BorderCase == -1;
    %         
    %         for n=[-1 0 1] 
    %             if round(x{j}(i,5))==2 || round(x{j}(i,5))==1
    %                 for k=[-1 0 1]
    %                     II{n+2,k+2,i}=ydatacrpd{i}(round(x{j}(i,4))+k,2+n); %2 indicates minimum distance
    %                 end                                                     %from left border                                    
    %             elseif round(x{j}(i,5))==3 || round(x{j}(i,5))==4
    %                 for k=[-2 -1 0 1 2]
    %                     II{n+2,k+3,i}=ydatacrpd{i}(round(x(i,4))+k,2+n); 
    %                 end
    %             end
    %         end
    % 
    %     %For BorderCase == 1
    %     elseif (round(x{j}(i,3))==2 || round(x{j}(i,3))==1) && BorderCase == 1;
    %         
    %         for n=[-1 0 1] 
    %             if round(x{j}(i,5))==2 || round(x{j}(i,5))==1
    %                 for k=[-1 0 1]
    %                     II{n+2,k+2,i}=ydatacrpd{i}(round(x{j}(i,4))+k,Size{i,j}(2)-2+n); 
    % 
    %                 end                                                                                         
    %             elseif round(x{j}(i,5))==3 || round(x{j}(i,5))==4
    %                 for k=[-2 -1 0 1 2]
    %                     II{n+2,k+3,i}=ydatacrpd{i}(round(x{j}(i,4))+k,Size{i,j}(2)-2+n); 
    %                 end
    %             end
    %         end
    %     
    %     % For sigmaX==3 or sigmaX==4 and Bordercase == 0
    %     elseif (round(x{j}(i,3))==3 || round(x{j}(i,3))==4) && BorderCaseWide == 0;
    %         
    %         for n=[-2 -1 0 1 2] 
    %             if round(x{j}(i,5))==2 || round(x{j}(i,5))==1
    %                 for k=[-1 0 1]
    %                     II{n+3,k+2,i}=ydatacrpd{i}(round(x{j}(i,4))+k,round(x{j}(i,2))+n);
    %                 end
    %             elseif round(x{j}(i,5))==3 || round(x{j}(i,5))==4
    %                 for k=[-2 -1 0 1 2]
    %                     II{n+3,k+3,i}=ydatacrpd{i}(round(x{j}(i,4))+k,round(x{j}(i,2))+n);
    %                 end
    %             end
    %         end
    %         
    %     elseif (round(x{j}(i,3))==3 || round(x{j}(i,3))==4) && BorderCaseWide == -1;
    %         
    %         for n=[-2 -1 0 1 2] 
    %             if round(x{j}(i,5))==2 || round(x{j}(i,5))==1
    %                 for k=[-1 0 1]
    %                     II{n+3,k+2,i}=ydatacrpd{i}(round(x{j}(i,4))+k,3+n);
    %                 end
    %             elseif round(x{j}(i,5))==3 || round(x{j}(i,5))==4
    %                 for k=[-2 -1 0 1 2]
    %                     II{n+3,k+3,i}=ydatacrpd{i}(round(x{j}(i,4))+k,3+n);
    %                 end
    %             end
    %         end
    %         
    %     elseif (round(x{j}(i,3))==3 || round(x{j}(i,3))==4) && BorderCaseWide == 1;
    %         
    %         for n=[-2 -1 0 1 2] 
    %             if round(x{j}(i,5))==2 || round(x{j}(i,5))==1
    %                 for k=[-1 0 1]
    %                     II{n+3,k+2,i}=ydatacrpd{i}(round(x{j}(i,4))+k,Size{i,j}(2)-3+n);
    %                 end
    %             elseif round(x{j}(i,5))==3 || round(x{j}(i,5))==4
    %                 for k=[-2 -1 0 1 2]
    %                     II{n+3,k+3,i}=ydatacrpd{i}(round(x{j}(i,4))+k,Size{i,j}(2)-3+n);
    %                 end
    %             end
    %         end
    %     end
    %     
    %         SII=[II{:,:,i}];
    %         x{j}(i,6)=sum(SII(:));
    %         XNorm{j}(i,6)=sum(SII(:));
    %     end
    % end



    %% plot

    % for i=1
    %     for j=1
    % 
    % figure(1)
    % hold on
    % imagesc(ydatacrpdR1{i,j})
    % plot(x{j}(i,2),x{j}(i,4),'r+')
    % hold off
    % axis([0 Size{i,j}(2)+1 0 size(ydatacrpd{i},1)+1])
    % 
    % xdatafit = linspace(-Size{i,j}(2),Size{i,j}(2)*2,300);
    % hdatafit = x{j}(i,1)*exp(-(xdatafit-x{j}(i,2)).^2/(2*x{j}(i,3)^2));
    % vdatafit = x{j}(i,1)*exp(-(xdatafit-x{j}(i,4)).^2/(2*x{j}(i,5)^2));
    % 
    % 
    % figure(2)
    % hold on
    % plot(ydatacrpdR1{i,j}(round(x{j}(i,4)),:),'ob');
    % plot(ydatacrpdR1{i,j}(:,round(x{j}(i,2))),'or');
    % plot(xdatafit,hdatafit,'-b',xdatafit,vdatafit,'-r')
    % axis([-5 30 0 1000])
    % hold off
    % legend('Horizontal','Vertical')
    %     end
    % end

    %% Save results
    %   save(strcat(Mainfolder,'DataMULTI/',num2str(Cell)),'x','XNorm','NSpots','SNR','ydatacrpd');

    % if ~exist(strcat(lionval.bacfolder,'Results'))
    %     mkdir(strcat(lionval.bacfolder,'Results'));
    % end

            display(strcat('Cell ',num2str(Cell),' analyzed'));
            display('Saving..');
        if ~exist(strcat(lionval.Mainfolder,'Results'),'dir')
            mkdir(strcat(lionval.Mainfolder,'Results'));
        end
        
        if exist(strcat(lionval.Mainfolder,'Results',lionval.OSslash,thisbacfolder,'.mat'))
            disp('This bac series has already been saved. Saving will be skipped')
        else
            save(strcat(lionval.Mainfolder,'Results',lionval.OSslash,thisbacfolder),'x','XNorm','NSpots','SNR','ydatacrpd','pixels','imgflipped');
            display('Save complete.');
        end


end

toc