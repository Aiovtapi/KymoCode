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

%% 
function DataStruct = GaussFitSimedit_Agarolysis(init,chan,Bacpics,Bacmask,cells,frames,IPTP)


%% Inputs 

lionval = Agar2lion(init,chan);
lionval.viewbac = 1:cells;

PSFSigma=1;
iterations=10;
fittype=4;

if exist(strcat(lionval.datapath,lionval.OSslash,'Results.mat'),'file')
    load(strcat(lionval.datapath,lionval.OSslash,'Results.mat'));
end

fprintf('\nTigerCut analyzed cell:')
    
for Cell=lionval.viewbac;
    
% % %     imgflip = 1;        % If you want to align the CFP images.
% % %     imgflipped = 0;
% % %     flipnumber = 0; flipdisplay = 0;
% % %     
% % %     while imgflip == 1; 

        clear XNorm x NSpots SNR ydatacrpd pixels

        ClipFactor=1;
        GaussFactor=1;

        thisbacfolder=strcat('Cell_',num2str(Cell,'%03.0f')); % Tif stack 
        bacseriepth=strcat(lionval.Mainfolder,lionval.OSslash,thisbacfolder,lionval.OSslash);
% 
%         
% % % %         % Check for whether flipping is needed based on the difchannel for
% % % %         % non difchannels. 
% % % %         if ~strcmp(lionval.viewchan,lionval.difchan)
% % % %             if exist(strcat(lionval.diffolder,'Results',lionval.OSslash,thisbacfolder,'.mat'))
% % % %                 flipimage=load(strcat(lionval.diffolder,'Results',lionval.OSslash,thisbacfolder,'.mat'),'imgflipped');
% % % %                 if flipimage.imgflipped==1
% % % %                     disp(['Flipping ',thisbacfolder])
% % % %                     GaussImRotate(bacseriepth);  
% % % %                 end
% % % %             else
% % % %                 disp('Dif-channel not computed. Cell orientation not defined. Do you want to continue? Click to continue')
% % % %                 waitforbuttonpress; close all
% % % %             end
% % % %             imgflip = 0;
% % % %             imgflipped = flipimage.imgflipped;           
% % % %         end
%         
% %         
%         d1{1}=readtimeseries(strcat(bacseriepth,'.tif'),'tif');
% 
%         data=dip_array(d1{1}); %turn into uint16 array
        bacsize = size(Bacpics{Cell,1});
        data = uint16(zeros(bacsize(1),bacsize(2),frames));
        
        for frami = 1:frames;
            data(:,:,frami) = Bacpics{Cell,frami};
        end

        %% Ze Defs

        Nspots=1; %ignore this name

        Tsize=frames;

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

        Tolerance=2;
        sigmachange=0.9;
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

            % Microfluidic channel boundaries (px)
            lowerboundchannel=6;
            higherboundchannel=18;
            
            % Intensity thresholding for outliers
            Data=ydatacrpd{i}; %%%(lowerboundchannel:higherboundchannel,:);
            [flag,A]=DetermineOutliers(Data,Tolerance,sigmachange,how,show);
            Outliersdata=~flag.*ydatacrpd{i}; %%%(lowerboundchannel:higherboundchannel,:);
            IntensityPeakThreshold = mean(nonzeros(Outliersdata(:)))+std(nonzeros(Outliersdata(:)));

            ydatacrpdR1{i,1}=Outliersdata;
            
            j=1;


            [x0{j}(i,:),Case{j}(i),ydatacrpdR1{i,j+1},Ydata{i,j},Size{i,j},Yg(i,j),Xg(i,j)]= ... 
                LionSpotter(ydatacrpdR1{i,j},SA,Sx,Sy,Px,Py,Bs,lob,upb);
            
           %still needs background level to be mediated from within the channel.   
            while x0{j}(i,1)>IPTP*IntensityPeakThreshold && j<20 && Size{i,j}(1)>0

            Data=mat2im(Ydata{i,j});
            
            % Main Fitting Algorithm
            
            [P,CRLB,LL,t]=gaussmlev2(Data,PSFSigma,iterations,fittype);
            
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
            
%             % define upper boundary (Size{i,j}(2) is the X size of the crpd-image)
%             ub = [realmax('double'),Size{i,j}(2),(Size{i,j}(2))^2,Size{i,j}(2), ...
%                 (Size{i,j}(2))^2]; %upper bounds

            % Do the Gaussian fitting

%             [x{j}(i,:),resnorm,residual,exitflag] = lsqcurvefit(@GaussPlosFunc, ...
%                x0{j}(i,:),xdata(:,j),Ydata{i,j},lb,ub,OPTIONS);


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
            
% % %            % Check for whether flipping is necessary. This will be skipped
% % %            % for non dif-channels, as flipping has already been done. 
% % %            if imgflip==1
% % %                 imgflip = x{1}(1,2)>size(ydatacrpd{1},2)/2;
% % %            end
% % %            
% % %            % Prevent flipping loops when the estimation of the spot is in
% % %            % the middle of the cell.
% % %            if flipnumber > 1;
% % %                 if flipdisplay == 0
% % %                     disp(strcat('cell',num2str(Cell),' is in a flipping loop, check bacseries in ImageJ.'))
% % %                     flipdisplay = 1;
% % %                 end
% % %                 imgflip = 0;
% % %            end
% % %            
% % %            if imgflip == 1
% % %                imgflipped = 1;
% % %                break
% % %            end



        %   Size{i,j}(1) is the height of the image.
        %   Size{i,j}(2) is the width of the image.

%             XNorm{j}(i,2)=x{j}(i,2)/Size{i,j}(2);
%             XNorm{j}(i,4)=(x{j}(i,4))/Size{i,j}(1);

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

% % %             % Breaking the loop for image flipping. 
% % %             if imgflip==1
% % %                 break
% % %             end
        end
        
% % %         % Flipping the images
% % %         if imgflip==1
% % %             disp(['Flipping ',thisbacfolder])
% % %             GaussImRotate(bacseriepth);
% % %             flipnumber = flipnumber + 1;
% % %         end
% % %     end

    NSpots=size(x,2);
%     LXNormOne=size(XNorm{1},1);
    LxOne=size(x{1},1);

    % make sure that each spot data is same length
    % because some spots will be there for just a few frames
    % and their last will be the length of the array.

    for j=2:NSpots;
        x{j}=[x{j};zeros(LxOne-size(x{j},1),5)];
%         XNorm{j}=[XNorm{j};zeros(LXNormOne-size(XNorm{j},1),4)];
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
         Cellidx = find(Bacmask{Cell,i});
         
         FCII{i}=ydatacrpd{i}; %%%(lowerboundchannel:higherboundchannel,:);
         
         FCIIoutput = sum(ydatacrpd{i}(Cellidx));

    end



    %% Save results
    
        % Display celli number
        if Cell>1
            for j=0:log10(Cell-1)
                fprintf('\b');
            end
        end
        fprintf(num2str(Cell))
        
%             display('Saving..');
%         if ~exist(strcat(lionval.Mainfolder,lionval.OSslash,'Results'),'dir')
%             mkdir(strcat(lionval.Mainfolder,lionval.OSslash,'Results'));
%         end
%         
%         if exist(strcat(lionval.Mainfolder,lionval.OSslash,'Results',lionval.OSslash,thisbacfolder,'.mat'),'file')
%             disp('This bac series has already been saved. Saving will be skipped')
%         else
%             save(strcat(lionval.Mainfolder,lionval.OSslash,'Results',lionval.OSslash,thisbacfolder),'x','NSpots','SNR','ydatacrpdR1','pixels'); %,'imgflipped');
%             display('Save complete.');
%         end
        
        DataStruct(chan,Cell).x = x;
        DataStruct(chan,Cell).SNR = SNR;
        DataStruct(chan,Cell).ydatacrpdR1 = ydatacrpdR1;
        DataStruct(chan,Cell).pixels = pixels;


end

save(strcat(lionval.datapath,lionval.OSslash,'Results.mat'),'DataStruct')

end