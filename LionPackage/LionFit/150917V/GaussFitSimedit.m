function GaussFitSimedit(user,exp,IPTPvalue)

%% Fit to spot 
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


%% 

tic

%% Inputs 
cropindx=1;
PSFSigma=1;
iterations=10;
fittype=4;


%%

initval=A001_Images_Set_Experiment(user,exp);

BacMainfolder = strcat(initval.basepath,initval.FiguresFolder,'BacPics',initval.OSslash);
Resultsfolder = strcat(initval.basepath,'Results');
if ~exist(Resultsfolder,'dir')
    mkdir(Resultsfolder)
end

FLchans = numel(initval.viewchan);
MDchans = initval.channelno;
finddif = strcmp(initval.difchan,initval.viewchan);
viewFLchans = {initval.viewchan{finddif},initval.viewchan{~finddif}};
difselected = any(finddif);
if nargin<3
    IPTPvalue = ones(1,FLchans);
end

if ~difselected
    choice = questdlg('Difchannel not selected. Rotation will not be calculated. Do you wish to continue?',...
        'Warning','Yes','No','No');
    switch choice
        case 'no'
            return
    end
end

%   FITTING PARAMETER INPUTS:
%   data:       PxPxN stack of images where N is number of PxP images
%   PSFSigma:   Microscope PSF sigma (sigma_0 for z fit, starting sigma for sigma fit)
%   iterations: Number of iterations (default=10)
%   fittype:    (default=1)
%    1: XY Position, Photons, Background
%    2: XY Position, Photons, Background, PSF Sigma
%    3: XYZ Position, Photons, Background 
%    4: XY Position, Photons, Background, PSF Sigma_X, PSF Sigma_Y


for mdchan = 1:MDchans
    disp(strcat('----- Multifluidic Device channel ',num2str(mdchan),' out of ',num2str(MDchans),' -----'))
    MDChanMainfolder=strcat(BacMainfolder,'Channel_0',num2str(mdchan),initval.OSslash);
    for flchan = 1:FLchans
        disp(strcat('----- LionFit on ',viewFLchans{flchan},' Channel -----'))
        FLChanMainfolder=strcat(MDChanMainfolder,viewFLchans{flchan},initval.OSslash);
        DifChanMainfolder=strcat(MDChanMainfolder,initval.difchan,initval.OSslash);
        
        thisdir = dir(FLChanMainfolder);
        thisdir = thisdir(arrayfun(@(x) ~strcmp(x.name(1),'.'),thisdir));
        Ncells = 1;
    
        for celli=1:Ncells

            imgflip = 1;        % If you want to align the CFP images.
            imgflipped = 0;
            flipnumber = 0; fliploop = 0;

            while imgflip == 1; 

                clear XNorm x NSpots SNR ydatacrpd pixels

                ClipFactor=1;
                GaussFactor=1;

                thisbacfolder=strcat('cell',num2str(celli,'%03.0f'));
                bacseriepth=strcat(FLChanMainfolder,thisbacfolder,initval.OSslash);

                % Check for whether flipping is needed based on the difchannel for
                % non difchannels. 
                if ~strcmp(viewFLchans{flchan},initval.difchan)
                    if exist(strcat(Resultsfolder,initval.OSslash,'Results_Ch',num2str(mdchan),'.mat'),'file')
                        load(strcat(Resultsfolder,initval.OSslash,'Results_Ch',num2str(mdchan),'.mat'));
                        imgflipped = DataStruct(1,celli).imgflipped;
                        if imgflipped==1
                            disp(['Flipping ',thisbacfolder])
                            GaussImRotate(bacseriepth);  
                        end
                    else
                        choice = questdlg('Difchannel not yet calculated. Rotation will not be applied. Do you wish to continue?',...
                            'Warning','Yes','No','No');
                        switch choice
                            case 'no'
                                return
                        end
                    end
                    imgflip = 0;
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

                Tolerance=3;
                sigmachange=0.2;
                how='positive';
                show=0;

                for i=1:Tsize

                    %Remove the x,y-Pads
                    ydata{i}=double(data(:,:,i));

                    %   Noticed that cropping causes shift in simulations!!
                    if cropindx==1;
                       [ydatacrpd{i},~]=Crop_Image(ydata{i});
                    else 
                        ydatacrpd{i}=ydata{i};
                    end
                    ydatacrpdR1{i,1}=ydatacrpd{i};

                    % Microfluidic channel boundaries (px)
                    lowerboundchannel=6;
                    higherboundchannel=18;

                    % Intensity thresholding for outliers
                    Data=ydatacrpdR1{i};
                    [flag,A]=DetermineOutliers(Data,Tolerance,sigmachange,how,show);
                    Outliersdata=~flag.*ydatacrpdR1{i};
                    IntensityPeakThreshold = mean(nonzeros(Outliersdata(:)))+std(nonzeros(Outliersdata(:)));

                    ydatacrpdR1{i,1}=Outliersdata;
                    

                    
                    j=1;


                    [x0{j}(i,:),Case{j}(i),ydatacrpdR1{i,j+1},Ydata{i,j},Size{i,j},Yg(i,j),Xg(i,j)]= ... 
                        LionSpotter(ydatacrpdR1{i,j},SA,Sx,Sy,Px,Py,Bs,lob,upb);

                   %still needs background level to be mediated from within the channel.   

                    while x0{j}(i,1)>IPTPvalue(flchan)*IntensityPeakThreshold && j<4 && Size{i,j}(1)>0

                        Fitdata=mat2im(Ydata{i,j});

                        % Main Fitting Algorithm

                        [P,CRLB,LL,t]=gaussmlev2(Fitdata,PSFSigma,iterations,fittype);

                        %       OUTPUTS:
                        %   P:      Found parameter values by fittype:
                        %    1: [X Y Photons BG] 
                        %    2: [X Y Photons BG Sigma] 
                        %    3: [X Y Photons BG Z]  
                        %    4: [X Y Photons BG SigmaX SigmaY] 
                        %   CRLB:   Cramer-Rao lower bound calculated using P
                        %   LL:     Log-Likelihood calculated using Stirling's approximation      
                        %   t:      Execution time of fitting code. 
                        %

                        P(1)=P(1)+1; %dip-image has origin (0,0)
                        P(2)=P(2)+1;

                        % transform to own format
                        x{j}(i,1)=P(3); 
                        x{j}(i,2)=P(1); 
                        x{j}(i,3)=P(5); 
                        x{j}(i,4)=P(2);
                        x{j}(i,5)=P(6);

                        % % define upper boundary (Size{i,j}(2) is the X size of the crpd-image)
                        % ub = [realmax('double'),Size{i,j}(2),(Size{i,j}(2))^2,Size{i,j}(2), ...
                        %    (Size{i,j}(2))^2]; %upper bounds
                        % 
                        % Do the Gaussian fitting
                        % 
                        % [x{j}(i,:),resnorm,residual,exitflag] = lsqcurvefit(@GaussPlosFunc, ...
                        %     x0{j}(i,:),xdata(:,j),Ydata{i,j},lb,ub,OPTIONS);
                        % 
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

                        % Prevent flipping loops when the estimation of the spot is in
                        % the middle of the cell.
                        if flipnumber > 1;
                           if fliploop == 0
                                disp(strcat('cell',num2str(celli),' is in a flipping loop, check bacseries in ImageJ.'))
                                fliploop = 1;
                           end
                           imgflip = 0;
                        end

                        if imgflip == 1
                            imgflipped = 1;
                            break
                        end

                        %   Size{i,j}(1) is the height of the image.
                        %   Size{i,j}(2) is the width of the image.

                        XNorm{j}(i,2)=x{j}(i,2)/Size{i,j}(2);

                        XNorm{j}(i,4)=x{j}(i,4)/Size{i,j}(1);

                        GaussmaskW=GaussFactor*(x{j}(i,3).^2+x{j}(i,5).^2)^(1/2);
                        ClipmaskR=ClipFactor*(x{j}(i,3).^2+x{j}(i,5).^2)^(1/2);

                        [~,~,ISPOT,Ibacklevel,spotim_clipped,bckim,ydatacrpdR1{i,j+1},pixels{j}(i)]=...
                            LionMasker(ydatacrpdR1{i,j},x{j}(i,2),x{j}(i,4),ClipmaskR,GaussmaskW);

                        if size(nonzeros(spotim_clipped(:)),1)<1
                           break
                        else
                           j=j+1;

                           [x0{j}(i,:),Case{j}(i),~,Ydata{i,j},Size{i,j},Yg(i,j),Xg(i,j)]= ... 
                           LionSpotter(ydatacrpdR1{i,j},SA,Sx,Sy,Px,Py,Bs,lob,upb);  
                        end

                    end
                    
                    if flipnumber > 1;
                        imgflip = 0;
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
                    flipnumber = flipnumber + 1;
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
            
            % Full Cell Intensity
            for i=1:Tsize;
                 FCII{i}=ydatacrpd{i}(lowerboundchannel:higherboundchannel,:);
                for j=1:NSpots
                    x{j}(i,7) = sum(sum(FCII{i}));
                end
            end

            %% plot

            for i=1
                for j=1
                    figure(1)
                    hold on
                    imagesc(ydatacrpd{i,j})
                    plot(x{j}(i,2),x{j}(i,4),'r+')
                    hold off
                    axis([0 Size{i,j}(2)+1 0 size(ydatacrpd{i},1)+1])
                end
            end

            %% Save results

            display(strcat('Cell ',num2str(celli),' analyzed'));
            display('Saving..');
%                 if ~exist(strcat(lionval.Mainfolder,'Results'),'dir')
%                     mkdir(strcat(lionval.Mainfolder,'Results'));
%                 end
% 
%                 if exist(strcat(lionval.Mainfolder,'Results',lionval.OSslash,thisbacfolder,'.mat'))
%                     disp('This bac series has already been saved. Saving will be skipped')
%                 else
%                     save(strcat(lionval.Mainfolder,'Results',lionval.OSslash,thisbacfolder),'x','XNorm','NSpots','ydatacrpd','pixels','imgflipped');
%                     display('Save complete in: ');
%                     display(strcat(lionval.Mainfolder,'Results',lionval.OSslash,thisbacfolder));
%                 end

            DataStruct(flchan,celli).x = x;
            DataStruct(flchan,celli).XNorm = XNorm;
            DataStruct(flchan,celli).NSpots = NSpots;
            DataStruct(flchan,celli).ydatacrpdR1 = ydatacrpdR1;
            DataStruct(flchan,celli).pixels = pixels;
            DataStruct(flchan,celli).imgflipped = imgflipped;
            DataStruct(flchan,celli).fliploop = fliploop;

            clear x XNorm NSpots ydatacrpdR1 pixels imgflipped fliploop

            save(strcat(Resultsfolder,initval.OSslash,'Results_Ch',num2str(mdchan),'.mat'),'DataStruct')
        end   
    end
    clear DataStruct
end
toc
end