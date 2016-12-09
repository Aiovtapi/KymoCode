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
function DataStruct = GaussFitSimedit_Agarolysis(init,chan,Bacpics,CBacmask,cells,frames,IPTP)


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
    
    clear XNorm x NSpots SNR ydatacrpd pixels

    ClipFactor=1;
    GaussFactor=2;

    bacsize = size(Bacpics{Cell,1});
    data = uint16(zeros(bacsize(1),bacsize(2),frames));

    for frami = 1:frames
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

    %% Load Data
    % parameters : [Amplitude, x0, sigmax, y0, sigmay, angel(in rad)]

    % Define some parameters 

    SA=3; % pixel step from maximum amplitude 
    Sx=3; Sy=3; % guess of sigmaX of spots
    Px=3; Py=3; % for optimal fitting add black 'pads' around detected spot
    Bs=1; Bss=2*Bs+1; % this defines size of black sqaure
    lob=1; upb=size(data(:,:,1),1); chanthickness=size(data(:,:,1),1); % define the boundaries of the channel

    % Determine outliers for determining intensity threshold.

    Tolerance=2;
    sigmachange=0.9;
    how='positive';
    show=0;

    for i=1:Tsize

        % Remove the x,y-Pads
        ydata{i}=double(data(:,:,i));

        % Noticed that cropping causes shift in simulations!!
        if lionval.cropindx==1
           [ydatacrpd{i},~]=Crop_Image(ydata{i});
        else 
            ydatacrpd{i}=ydata{i};
        end
        ydatacrpdR1{i,1}=ydatacrpd{i};

        % Intensity thresholding for outliers
        Data=ydatacrpd{i};
        [flag,A]=DetermineOutliers(Data,Tolerance,sigmachange,how,show);
        Outliersdata=~flag.*ydatacrpd{i};
        IntensityPeakThreshold = mean(nonzeros(Outliersdata(:)))+std(nonzeros(Outliersdata(:)));

        ydatacrpdR1{i,1}=Outliersdata;

        j=1;

        [x0{j}(i,:),Case{j}(i),ydatacrpdR1{i,j+1},Ydata{i,j},Size{i,j},Yg(i,j),Xg(i,j)]= ... 
            LionSpotter(ydatacrpdR1{i,j},SA,Sx,Sy,Px,Py,Bs,lob,upb);

       % still needs background level to be mediated from within the channel.   
        while x0{j}(i,1)>IPTP*IntensityPeakThreshold && j<20 && Size{i,j}(1)>0

            Data=mat2im(Ydata{i,j});

            % Main Fitting Algorithm
            [P,CRLB,LL,t]=gaussmlev2(Data,PSFSigma,iterations,fittype);
            
            % OUTPUTS:
            %   P:      Found parameter values by fittype:
            %    1: [X Y Photons BG] 
            %    2: [X Y Photons BG Sigma] 
            %    3: [X Y Photons BG Z]  
            %    4: [X Y Photons BG SigmaX SigmaY] 
            %   CRLB:   Cramer-Rao lower bound calculated using P
            %   LL:     Log-Likelihood calculated using Stirling's approximation      
            %   t:      Execution time of fitting code. 
            
            P(1)=P(1)+1; %dip-image has origin (0,0)
            P(2)=P(2)+1;
            
            % transform to own format
            x{j}(i,1)=P(3); 
            x{j}(i,2)=P(1); 
            x{j}(i,3)=P(5); 
            x{j}(i,4)=P(2);
            x{j}(i,5)=P(6);
          
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
    end
        

    NSpots=size(x,2);
%   LXNormOne=size(XNorm{1},1);
    LxOne=size(x{1},1);

    % make sure that each spot data is same length
    % because some spots will be there for just a few frames
    % and their last will be the length of the array.

    for j=2:NSpots
        x{j}=[x{j};zeros(LxOne-size(x{j},1),5)];
%       XNorm{j}=[XNorm{j};zeros(LXNormOne-size(XNorm{j},1),4)];
    end


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
         Cellidx = find(CBacmask{Cell,i});
         FCII{i}=ydatacrpd{i}; 

        for j=1:NSpots
            %This is the full cell integrated intensity
            x{j}(i,7) = sum(FCII{i}(Cellidx));
            
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

%             ROI{i,j}=interp2(XX,YY,ydatacrpdR1{i,j},xx0,yy0,'cubic');
                %Have to fix the unsupport for interp2
                
              ROI{i,j}=ydatacrpdR1{i,j};
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



    %% Save results
    
        % Display celli number
        if Cell>1
            for j=0:log10(Cell-1)
                fprintf('\b');
            end
        end
        fprintf(num2str(Cell))
        
        DataStruct(chan,Cell).x = x;
        DataStruct(chan,Cell).SNR = SNR;
        DataStruct(chan,Cell).ydatacrpdR1 = ydatacrpdR1;
        DataStruct(chan,Cell).pixels = pixels;


end

save(strcat(lionval.datapath,lionval.OSslash,'Results.mat'),'DataStruct')

end