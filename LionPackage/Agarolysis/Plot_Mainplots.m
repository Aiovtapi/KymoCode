clear 
% close all
clc

%%
folder = '/Users/rleeuw/Work/Data/141230_Kymo_Results';
slash = '/';
screensize = get( groot, 'Screensize' );
screensize(4) = screensize(4);
screensize2 = screensize;
screensize2(3) = screensize(3)/8*3;
screensize3 = screensize;
screensize3(3) = screensize(3)/8*5;
pos = [225   327   823   467];
Font = 'Raleway Medium';
AW = 1.5;       % axis thickness
FS = 14;    FS2 = 18;       % font size
LW = 1.2;       % Line width
green = [0, 0.5, 0];
blue = [0, 0.5, 1];
red = [0.8,0,0];
oran = [1, 0.5, 0];

Title = {'A.','B.'};
TitleFont = 'RalewaySemiBold';
TitleFS = 20;   TitleFS2 = 26;
Tx = -0.1; Ty = 0.05;

px=0.1589;
kdnancal = 1100;
ktuscal = 600;
adnancal = 990;
atuscal = 3300;

CLbound = [14 39];

CL = 'Cell Length (\mum)';
CC = 'Replication Cycle';
P = 'Normalized Position Within Cell';
SI = 'Spot Intensity';
SIN = 'Spot Intensity / Number of Spots';
NS = 'Number of Spots';



%%
k_getdata = true;
a_getdata = true;


k_PvsCL = false;
k_IvsP = false;
k_NSvsP = false;
k_Heatmap = false;
k_FCIvsCL = false;
k_royplots = false;
k_numspots =false;
k_numprots = false;
k_numprots2 = true;

a_PvsCL = false;
a_IvsP = false;
a_NSvsP = false;
a_Heatmap = false;
a_FCIvsCL = false;
a_royplots = false;
a_numspots = false;
a_numprots = false;

c_PvsCL = false;
c_IvsP = false;
c_NSvsP = false;
c_Heatmap = false;

y_PvsCL = false;
y_IvsP = false;
y_NSvsP = false;
y_Heatmap = false;

r_PvsCL = false;
r_IvsP = false;
r_NSvsP = false;
r_Heatmap = false;

all_Heatmap = false;

% c_PvsCL = true;
% c_IvsP = true;
% c_NSvsP = true;
% c_Heatmap = true;
% 
% y_PvsCL = true;
% y_IvsP = true;
% y_NSvsP = true;
% y_Heatmap = true;
% 
% r_PvsCL = true;
% r_IvsP = true;
% r_NSvsP = true;
% r_Heatmap = true;


%% Kymo
    % Get Data 
    if k_getdata == 1

        k_channels = [1,2];
        k_Intensityval = [60,60,60];

        kLcfp=[];    kLyfp=[];    kLrfp=[];
        kPcfp=[];    kPyfp=[];    kPrfp=[];
        kIcfp=[];    kIyfp=[];    kIrfp=[];
        kFcfp=[];    kFyfp=[];    kFrfp=[];
        kncfp=[];    knyfp=[];    knrfp=[];
        kNSc=[];     kNSy=[];     kNSr=[];
        kCLc=[];     kCLy=[];     kCLr=[];
        kIc=[];      kIy=[];      kIr=[];
        kFc=[];      kFy=[];      kFr=[];
        kQc=[];      kQy=[];      kQr=[];
        
        allframes = [];
        
        YFP.filterval=700;
        RFP.filterval=700;
        CFP.filterval=1000;

        for i=k_channels;
            E{i}=load(strcat(folder,slash,'Results_Ch',num2str(i),'.mat')); 
        end

        allCFP_L = [];
        % i: multiluidic channel
        % j: cell number
        % h: frame number
        % k: spot number
        
        for i=k_channels

            Ncells{i}=size(E{i}.DataStruct,2);


            for j=1:Ncells{i}     

                CFPx{i}{j}=E{i}.DataStruct(1,j).x;
                YFPx{i}{j}=E{i}.DataStruct(2,j).x;
                RFPx{i}{j}=E{i}.DataStruct(3,j).x;

                NspotsCFP=size(CFPx{i}{j},2);
                NspotsYFP=size(YFPx{i}{j},2);
                NspotsRFP=size(RFPx{i}{j},2);

                frames{i}{j} = size(CFPx{i}{j}{1},1);
                allframes = [allframes,frames{i}{j}];
                
                if NspotsCFP==0
                    CFPx{i}{j}{1}=[];
                else
                    for h = 1:frames{i}{j}
                        Ic = 0;
                        for k=1:NspotsCFP
                            length = size(E{i}.DataStruct(1,j).ydatacrpdR1{h,k},2);
                            kLcfp=[kLcfp h/frames{i}{j}];
                            kPcfp=[kPcfp CFPx{i}{j}{k}(h,2)/length];
                            kIcfp=[kIcfp CFPx{i}{j}{k}(h,1)];
                            IcfpM{i}{j}(h,k)=CFPx{i}{j}{k}(h,1);
                            IcfpM{i}{j}(IcfpM{i}{j}<CFP.filterval)=0;
                            kFcfp=[kFcfp CFPx{i}{j}{k}(h,7)];
                            Ic = Ic + CFPx{i}{j}{k}(h,1);
                        end
                        thisy = E{i}.DataStruct(1,j).ydatacrpdR1(h,:);
                        numspot = sum(~cellfun('isempty',thisy));
%                         clength = size(thisy{1},2);
                        clength = h/frames{i}{j};
                        kNSc = [kNSc,numspot];
                        kCLc = [kCLc,clength];
                        kIc = [kIc, Ic];
                        kFc = [kFc, CFPx{i}{j}{1}(h,7)];
                        kQc = [kQc, Ic/CFPx{i}{j}{1}(h,7)*100];
                    end
                end


                if NspotsYFP==0
                    YFPx{i}{j}{1}=[];
                else
                    for h = 1:frames{i}{j}
                        Iy = 0;
                        for k=1:NspotsYFP
                            length = size(E{i}.DataStruct(2,j).ydatacrpdR1{h,k},2);
                            kLyfp=[kLyfp h/frames{i}{j}];
                            lyfp{i}{j}(h)=h/frames{i}{j};
                            kPyfp=[kPyfp YFPx{i}{j}{k}(h,2)/length];
                            kIyfp=[kIyfp YFPx{i}{j}{k}(h,1)];
                            IyfpM{i}{j}(h,k)=YFPx{i}{j}{k}(h,1);
                            IyfpM{i}{j}(IyfpM{i}{j}<YFP.filterval)=0;
                            kFyfp=[kFyfp YFPx{i}{j}{k}(h,7)];
                            Iy = Iy + YFPx{i}{j}{k}(h,1);
                        end
                        thisy = E{i}.DataStruct(2,j).ydatacrpdR1(h,:);
                        numspot = sum(~cellfun('isempty',thisy));
%                         clength = size(thisy{1},2);
                        clength = h/frames{i}{j};
                        kNSy = [kNSy,numspot];
                        kCLy = [kCLy,clength];
                        kIy = [kIy, Iy];
                        kFy = [kFy, YFPx{i}{j}{1}(h,7)];
                        kQy = [kQy, Iy/YFPx{i}{j}{1}(h,7)*100];
                    end
                end


                if NspotsRFP==0
                    RFPx{i}{j}{1}=[];
                else
                    for h = 1:frames{i}{j}
                        Ir = 0;
                        for k=1:NspotsRFP
                            length = size(E{i}.DataStruct(3,j).ydatacrpdR1{h,k},2);
                            kLrfp=[kLrfp h/frames{i}{j}];
                            lrfp{i}{j}=h/frames{i}{j};
                            kPrfp=[kPrfp RFPx{i}{j}{k}(h,2)/length];
                            kIrfp=[kIrfp RFPx{i}{j}{k}(h,1)];                   
                            IrfpM{i}{j}(h,k)=RFPx{i}{j}{k}(h,1);
                            IrfpM{i}{j}(IrfpM{i}{j}<RFP.filterval)=0;
                            kFrfp=[kFrfp RFPx{i}{j}{k}(h,7)];
                            Ir = Ir + RFPx{i}{j}{k}(h,1);
                        end
                        thisy = E{i}.DataStruct(3,j).ydatacrpdR1(h,:);
                        numspot = sum(~cellfun('isempty',thisy));
%                         clength = size(thisy{1},2);
                        clength = h/frames{i}{j};
                        kNSr = [kNSr,numspot];
                        kCLr = [kCLr,clength];
                        kIr = [kIr, Ir];
                        kFr = [kFr, RFPx{i}{j}{1}(h,7)];
                        kQr = [kQr, Ir/RFPx{i}{j}{1}(h,7)*100];
                    end
                end

            end
        end

        Yremove = kIyfp == 0;
        kIyfp(Yremove) = [];
        kLyfp(Yremove) = [];
        kPyfp(Yremove) = [];
        kFyfp(Yremove) = [];

        Rremove = kIrfp == 0;
        kIrfp(Rremove) = [];
        kLrfp(Rremove) = [];
        kPrfp(Rremove) = [];
        kFrfp(Rremove) = [];

        Cremove = kIcfp == 0;
        kIcfp(Cremove) = [];
        kLcfp(Cremove) = [];
        kPcfp(Cremove) = [];
        kFcfp(Cremove) = [];
        
        YR2 = kIy == 0;
        kNSy(YR2) = [];
        kCLy(YR2) = [];
        kIy(YR2)  = [];
        kFy(YR2)  = [];
        kQy(YR2)  = [];

        RR2 = kIr == 0;
        kNSr(RR2) = [];
        kCLr(RR2) = [];
        kIr(RR2)  = [];
        kFr(RR2)  = [];
        kQr(RR2)  = [];

        CR2 = kIc == 0;
        kNSc(CR2) = [];
        kCLc(CR2) = [];
        kIc(CR2)  = [];
        kFc(CR2)  = [];
        kQc(CR2)  = [];

        %Total spot intensities
        %Spot number

        for i=k_channels
            for j=1:Ncells{i}

                IyfpSpotTotal{i,j}=sum(IyfpM{i}{j},2);
                IcfpSpotTotal{i,j}=sum(IcfpM{i}{j},2);
                IrfpSpotTotal{i,j}=sum(IrfpM{i}{j},2);

                %find number of spots
                %find nonzero elements
                [RowYfpnz,ColYfpnz]=find(IyfpM{i}{j}>0);
                %find multitude of spots by sorting and counting occurences by
                %sorting
                [MultiSpotYFP,SortIdxYfp]=sort(RowYfpnz);         
                %sort, occurrence, then unique and link
                DD=MultiSpotYFP;
                [aanz,bbnz]=hist(DD,unique(DD)); % aa is multitude, bb is row number

                % now we now each row's multitude, thus per frame the number of spots. 
                % Per cell j we need to create one array with spot numbers from
                % first to second frame. So we need to fill in the gaps because we
                % only have the nonzero elements now 
                YFPSpotNumber{i,j}=zeros(1,frames{i}{j});
                for L=1:size(bbnz,1)
                YFPSpotNumber{i,j}(bbnz(L))=aanz(L);
                end

                % do same for RFP
                [RowRfpnz,ColRfpnz]=find(IrfpM{i}{j}>0);
                [MultiSpotRFP,SortIdxRfp]=sort(RowRfpnz);         
                DD=MultiSpotRFP;
                [aanz,bbnz]=hist(DD,unique(DD)); 

                RFPSpotNumber{i,j}=zeros(1,frames{i}{j});

                for L=1:size(bbnz,1)
                RFPSpotNumber{i,j}(bbnz(L))=aanz(L);
                end
            end
        end
    end

    % Position vs. cell length
    if k_PvsCL==1

        figure
        set(gcf, 'position', screensize)
        set(gcf,'Color',[0.95 0.95 0.95])

        % CFP
        subplot(1,3,1)
        hold on
        scatter(single(kLcfp),kPcfp,kIcfp,blue,'filled');
        % myfit=polyfit(Acfp,Bcfp,4);
        % x=0:0.01:1;
        % y=polyval(myfit,x);
        % plot(x,y,'r','LineWidth',5)
        hold off
        ylabel(P); 
        t=title('Kymo data: CFP');
        axis([0 1 0 1])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font,'FontSize',FS)

        % YFP
        subplot(1,3,2)
        hold on
        scatter(single(kLyfp),kPyfp,kIyfp,green,'filled');
        % myfit=polyfit(Ayfp,Byfp,4);
        % x=0:0.01:1;
        % y=polyval(myfit,x);
        % plot(x,y,'k','LineWidth',5)
        hold off
        xlabel(CL);
        t=title('Kymo data: YFP');
        axis([0 1 0 1])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font,'FontSize',FS)

        % RFP
        subplot(1,3,3)
        hold on
        scatter(single(kLrfp),kPrfp,kIrfp,'r','filled');
        % myfit=polyfit(Arfp,Brfp,4);
        % x=0:0.01:1;
        % y=polyval(myfit,x);
        % plot(x,y,'k','LineWidth',5)
        hold off
        t=title('Kymo data: RFP');
        axis([0 1 0 1])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font,'FontSize',FS)
    end

    % Intensity vs. position
    if k_IvsP ==1

        figure
        set(gcf, 'position', screensize)
        set(gcf,'Color',[0.95 0.95 0.95])

        % CFP
        subplot(1,3,1)
        hold on
        scatter(kPcfp,kIcfp,'x','MarkerEdgeColor',blue);
        myfit=polyfit(kPcfp,kIcfp,4);
        x=0:00.1:1;
        y=polyval(myfit,x);
        plot(x,y,'k','LineWidth',3)
        ylabel(SI); 
        t=title('Kymo data: CFP');
        hold off
        axis([0 1 -0.1 35])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font,'FontSize',FS)

        % YFP
        subplot(1,3,2)
        hold on
        scatter(kPyfp,kIyfp,'x','MarkerEdgeColor',green);
        myfit=polyfit(kPyfp,kIyfp,4);
        x=0:00.1:1;
        y=polyval(myfit,x);
        plot(x,y,'k','LineWidth',3)
        xlabel('Position in cell');
        t=title('Kymo data: YFP');
        hold off
        axis([0 1 -0.1 35])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font,'FontSize',FS)

        % RFP
        subplot(1,3,3)
        hold on
        scatter(kPrfp,kIrfp,'r','x');
        myfit=polyfit(kPrfp,kIrfp,4);
        x=0:00.1:1;
        y=polyval(myfit,x);
        plot(x,y,'k','LineWidth',3)
        t=title('Kymo data: RFP');
        hold off
        axis([0 1 -0.1 35])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font,'FontSize',FS)
    end

    % Numspots vs. position
    if k_NSvsP==1

        figure
        set(gcf, 'position', screensize)
        set(gcf,'Color',[0.95 0.95 0.95])

        % CFP
        subplot(1,3,1)
        [numbin,edges] = histcounts(kPcfp,20);
        norm = max(numbin)/35;
        X = diff(edges);
        X = cumsum(X) - X(1)/2;
        hold on
        scatter(kPcfp,kIcfp,'x','MarkerEdgeColor',blue);
        plot(X,numbin/norm,'k','LineWidth',3)
        hold off
        ylabel(SIN);
        t=title('Kymo data: CFP');
        axis([0 1 -0.1 40])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font,'FontSize',FS)

        % YFP
        subplot(1,3,2)
        [numbin,edges] = histcounts(kPyfp,20);
        norm = max(numbin)/35;
        X = diff(edges);
        X = cumsum(X) - X(1)/2;
        hold on
        scatter(kPyfp,kIyfp,'x','MarkerEdgeColor',green);
        plot(X,numbin/norm,'k','LineWidth',3)
        hold off
        xlabel(P);
        t=title('Kymo data: YFP');
        axis([0 1 -0.1 40])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font,'FontSize',FS)

        % RFP
        subplot(1,3,3)
        [numbin,edges] = histcounts(kPrfp,20);
        norm = max(numbin)/35;
        X = diff(edges);
        X = cumsum(X) - X(1)/2;
        hold on
        scatter(kPrfp,kIrfp,'r','x');
        plot(X,numbin/norm,'k','LineWidth',3)
        hold off
        t=title('Kymo data: RFP');
        axis([0 1 -0.1 40])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font,'FontSize',FS)
    end

    % Heatmaps
    if k_Heatmap==1

        bins = 15;
        thisedge2{1} = linspace(min(kLcfp),max(kLcfp),bins+5);
        thisedge2{2} = (0:bins)/bins;

        % CFP
        kNumcfp(1,:) = kLcfp;
        kNumcfp(2,:) = kPcfp;
        
        [bin_cfp,idx_cfp]=find(kIcfp>3000);

        kNumCFP(:,idx_cfp)=kNumcfp(:,idx_cfp);
        kNumCFPnz(1,:)=nonzeros(kNumCFP(1,:));
        kNumCFPnz(2,:)=nonzeros(kNumCFP(2,:));

        figure(21)
        subplot(3,2,1)
        Heatmap = hist3(kNumCFPnz','Edges',thisedge2);
        pcolor(thisedge2{1},(thisedge2{2}),Heatmap');
        colormap(winter) % heat map
        t=title('Kymo data: CFP');
        grid on
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font,'FontSize',FS)

%         figure(22)
%         subplot(1,3,1)
%         hold on
%         hist3(kNumcfp','Edges',thisedge2)
%         colormap(hot) % heat map
%         set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
%         xlabel(CL); ylabel(P);zlabel(NS)
%         t=title('Kymo data: CFP');
%         grid off
%         hold off
%         axis([min(kLcfp),max(kLcfp),0,1])
%         view(3)
%         set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
%         set(t,'Fontname',Font,'FontSize',FS)

        % YFP
        kNumyfp(1,:) = kLyfp;
        kNumyfp(2,:) = kPyfp;

        %Filter on Intensity
        [bin_yfp,idx_yfp]=find(kIyfp>1500); %980=Noise-Fraction 0.57

        kNumYFP(:,idx_yfp)=kNumyfp(:,idx_yfp);
        kNumYFPnz(1,:)=nonzeros(kNumYFP(1,:));
        kNumYFPnz(2,:)=nonzeros(kNumYFP(2,:));

        FilteredSignalYFP=size(kNumYFPnz,2)/size(kIyfp,2); % Percentage of total signal

        figure(21)
        subplot(3,2,3)
        Heatmap = hist3(kNumyfp','Edges',thisedge2);
        pcolor(thisedge2{1},(thisedge2{2}),Heatmap');
        colormap(autumn) % heat map
        ylabel(P);
        t=title('Kymo data: YFP');
        grid on
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font,'FontSize',FS)

%         figure(22)
%         subplot(1,3,2)
%         hold on
%         hist3(kNumyfp','Edges',thisedge2)
%         colormap(hot) % heat map
%         set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
%         xlabel(CL); ylabel(P);
%         t=title('Kymo data: YFP');
%         grid off
%         hold off
%         axis([min(kLyfp),max(kLyfp),0,1])
%         view(3)
%         set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
%         set(t,'Fontname',Font,'FontSize',FS)

        % RFP
        kNumrfp(1,:) = kLrfp;
        kNumrfp(2,:) = kPrfp;
        
        [bin_rfp,idx_rfp]=find(kIrfp>2000); %980=Noise-Fraction 0.57

        kNumRFP(:,idx_rfp)=kNumrfp(:,idx_rfp);
        kNumRFPnz(1,:)=nonzeros(kNumRFP(1,:));
        kNumRFPnz(2,:)=nonzeros(kNumRFP(2,:));

        figure(21)
        subplot(3,2,5)
        Heatmap = hist3(kNumRFPnz','Edges',thisedge2);
        h = pcolor(thisedge2{1},(thisedge2{2}),Heatmap');
        colormap(hot) % heat map
        xlabel(CL);
        t=title('Kymo data: RFP');
        grid on
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font,'FontSize',FS)

%         figure(22)
%         subplot(1,3,3)
%         hold on
%         hist3(kNumrfp','Edges',thisedge2)
%         colormap(hot) % heat map
%         set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
%         xlabel(CL); ylabel(P);
%         t=title('Kymo data: RFP');
%         grid off
%         hold off
%         axis([min(kLrfp),max(kLrfp),0,1])
%         view(3)
%         set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
%         set(t,'Fontname',Font,'FontSize',FS)
    end

    % Full cell intensity vs. celllength
    if k_FCIvsCL==1
        clear plotcfp plotyfp plotrfp
        fig6 = figure(6);
        set(fig6,'Position',[20,300,1800,500])
        YFPcal=953;
        RFPcal=990;

        % YFP

        plotcfp(1,:) = kLcfp;
        plotcfp(2,:) = kFcfp;
        plotcfp(3,:) = kIcfp*k_Intensityval(1);
        plotcfp = unique(plotcfp','rows')';

        subplot(1,3,1)
        hold on
        scatter(plotcfp(1,:),plotcfp(2,:),'b','x');
        scatter(plotcfp(1,:),plotcfp(3,:),'r','x');
        myfit=polyfit(plotcfp(1,:),plotcfp(2,:),1);
        x=0:0.01:1;
        y=polyval(myfit,x);
        plot(x,y,'k','LineWidth',3)
        xlabel('Cell Length'); ylabel('Normalized full cell intensity'); 
        title('Kymo data: CFP')
        hold off
        % axis([0 1 -0.1 1.1])

        % YFP

        plotyfp(1,:) = kLyfp;
        plotyfp(2,:) = kFyfp/YFPcal;
        plotyfp(3,:) = kIyfp*k_Intensityval(2)/YFPcal;
        plotyfp = unique(plotyfp','rows')';

        subplot(1,3,2)
        hold on
        scatter(plotyfp(1,:),plotyfp(2,:),'b','x');
        % scatter(plotyfp(1,:),plotyfp(3,:),'r','x');
        for i=k_channels
            for j=1:Ncells{i}
        scatter((1:size(IyfpSpotTotal{i,j},1))/size(IyfpSpotTotal{i,j},1),IyfpSpotTotal{i,j}./YFP.filterval,'r')
            end
        end
        myfit=polyfit(plotyfp(1,:),plotyfp(2,:),1);
        x=0:0.01:1;
        y=polyval(myfit,x);
        plot(x,y,'k','LineWidth',3)
        xlabel('Normalized Cell Time (-)'); ylabel('Number of Tus (-)'); 
        title('Tus Stoichiometry vs. Time')
        hold off
        axis([0 1 0.01 30])
        set(gca,'FontSize',18)


        % RFP

        plotrfp(1,:) = kLrfp;
        plotrfp(2,:) = kFrfp/RFPcal;
        plotrfp(3,:) = kIrfp*k_Intensityval(3)/RFPcal;
        plotrfp = unique(plotrfp','rows')';

        subplot(1,3,3)
        hold on
        % scatter(plotrfp(1,:),plotrfp(2,:),'b','x');
        scatter(plotrfp(1,:),plotrfp(3,:),'r','x');
        myfit=polyfit(plotrfp(1,:),plotrfp(2,:),1);
        x=0:0.01:1;
        y=polyval(myfit,x);
        % plot(x,y,'k','LineWidth',3)
        xlabel('Cell Length'); ylabel('Normalized full cell intensity'); 
        title('Kymo data: RFP')
        hold off
        % axis([0 1 -0.1 1.1])
    end

    % YFP Spot number vs time
    if k_royplots==1
        figure

        YfpSpot=[];
        Frame=[];

        for i=k_channels
            for j=1:Ncells{i}
                YfpSpot=[YfpSpot YFPSpotNumber{i,j}];
                Frame=[Frame (1:frames{i}{j})/frames{i}{j}];
            end
        end

        for i=1:max(YfpSpot);
        p{i} = find(YfpSpot == i);
        m(i) = mean(Frame(p{i}));
        Spotstd(i)=std(Frame(p{i}));
        end

        Y=linspace(0,max(YfpSpot),max(YfpSpot));

        Xl=m-Spotstd;
        Xr=m+Spotstd;

        hold on
        scatter(Frame,YfpSpot,'r','filled')
        plot(m,Y,'b--','LineWidth',3)
        plot(Xl,Y,'b','LineWidth',1)
        plot(Xr,Y,'b','LineWidth',1)
        hold off
        title('Tus Spot number vs. Normalized Cell Cycle')
        xlabel('Normalized Cell Time (-)'); ylabel('Number of Spots (-)')
        set(gca,'FontSize',18)


        %%
        figure
        hold on

        for i=k_channels
            for j=1:Ncells{i}
                scatter((1:size(IrfpSpotTotal{i,j},1))/size(IrfpSpotTotal{i,j},1),IrfpSpotTotal{i,j}/RFP.filterval,'r','filled')
            end
        end
    end
    
    % Numspots
    if k_numspots==1
        minv=floor(min([min(kCLc),min(kCLy),min(kCLr)]));
        maxv=ceil(max([max(kCLc),max(kCLy),max(kCLr)]));
        binsize = 0.075;

        bins = minv:binsize:maxv;
        nbins = numel(bins)-1;
        xbin = bins(2:end)-binsize/2;

        [MeanMy, MeanMr, MeanMc, NumMy, NumMr, NumMc, StdEMy, StdEMr, StdEMc] = deal(zeros(nbins,1));

        for j = 1:nbins
            yindex = find(kCLy>bins(j)&kCLy<bins(j+1));
            rindex = find(kCLr>bins(j)&kCLr<bins(j+1));
            cindex = find(kCLc>bins(j)&kCLc<bins(j+1));

            thismy = mean(kNSy(yindex));
            thismr = mean(kNSr(rindex));
            thismc = mean(kNSc(cindex));

            MeanMy(j) = thismy;
            MeanMr(j) = thismr;
            MeanMc(j) = thismc;

            StdEMy(j) = std(kNSy(yindex)); %/sqrt(numel(yindex));
            StdEMr(j) = std(kNSr(rindex)); %sqrt(numel(rindex));
            StdEMc(j) = std(kNSc(cindex)); %sqrt(numel(cindex));

            NumMy(j) = numel(yindex);
            NumMr(j) = numel(rindex);
            NumMc(j) = numel(cindex);
        end
        nany = isnan(MeanMy);
        MeanMy(nany)=[]; NumMy(nany)=[]; StdEMy(nany)=[]; ybin = xbin(~nany);
        nanr = isnan(MeanMr);
        MeanMr(nanr)=[]; NumMr(nanr)=[]; StdEMr(nanr)=[]; rbin = xbin(~nanr);
        nanc = isnan(MeanMc);
        MeanMc(nanc)=[]; NumMc(nanc)=[]; StdEMc(nanc)=[]; cbin = xbin(~nanc);

        ss=0.01;

        figure
        hold on
        
        scatter(kCLy,kNSy,30,green,'x','MarkerEdgeAlpha',0.15,'MarkerFaceAlpha',0.15);
        errorbar(ybin,MeanMy,StdEMy,'vertical','LineStyle','none','LineWidth',LW,'Color',green)
        
%         scatter(kCLc,kNSc,30,blue,'x','MarkerEdgeAlpha',0.15,'MarkerFaceAlpha',0.15);
%         errorbar(cbin-ss,MeanMc,StdEMc,'vertical','LineStyle','none','LineWidth',LW,'Color',blue)
%         scatter(kCLr,kNSr,30,'rx','MarkerEdgeAlpha',0.15,'MarkerFaceAlpha',0.15);
%         errorbar(rbin+ss,MeanMr,StdEMr,'vertical','LineStyle','none','LineWidth',LW,'Color','r')

        fity = fit(ybin(1:end-3)',MeanMy(1:end-3),'poly1');
        fitplot = plot(fity,'g');

        legend off
%         axis([10*px 45*px 0 6])
        axis([0 1 0 6])

        scc = 2;
        h2=scatter(xbin,MeanMy,NumMy*scc,green,'LineWidth',LW);
        scatter(xbin,MeanMy,NumMy*scc,green,'filled','LineWidth',LW,'MarkerFaceAlpha',0.2)
        
%         h1=scatter(cbin-ss,MeanMc,NumMc*scc,blue,'LineWidth',LW);
%         scatter(cbin-ss,MeanMc,NumMc*scc,blue,'filled','LineWidth',LW,'MarkerFaceAlpha',0.2)
%         h2=scatter(rbin+ss,MeanMr,NumMr*scc,'r','LineWidth',LW);
%         scatter(rbin+ss,MeanMr,NumMr*scc,'r','filled','LineWidth',LW,'MarkerFaceAlpha',0.2)
        
        vlplot = vline(bins,'k:');
        set(vlplot,'LineWidth',0.5)
        vl = vline(0.9,'k--');
        set(vl,'LineWidth',LW)
        
        hold off
        xlabel(CC)
        ylabel('Number of spots')
        set(fitplot,'LineWidth',2)
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS,'Color',[0.95 0.95 0.95])
        set(gcf,'Position',pos)
        set(gca,'YGrid','on')
        set(gca,'GridLineStyle','-')
                legend(h2,'Tus')
                
        if true
            figure
            hold on
            X2=[ybin,fliplr(ybin)];
            Y1=[(MeanMy-StdEMy)',fliplr((MeanMy+StdEMy)')];
            fila = fill(X2,Y1,green);
            set(fila,'LineStyle','none','FaceAlpha','0.1')
            plot(ybin,MeanMy+StdEMy,'-','LineWidth',LW,'Color',green);
            plot(ybin,MeanMy-StdEMy,'-','LineWidth',LW,'Color',green);
            plot(ybin,MeanMy,'-','LineWidth',3*LW,'Color',green)        

            axis([0 1 -2 7])
            vl = vline(0.9,'k--');
            set(vl,'LineWidth',LW)
            hold off
            legend off
            xlabel(CC)
            ylabel('Proteins per cell')
            set(gcf,'Position',pos)
            set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS,'YGrid','on','GridLineStyle','-','Color',[0.95 0.95 0.95])
        end
                
%         legend([h1,h2,h3],{'dif','Tus','DnaN'})
    end
    
   % Numprots
    if k_numprots==1
        minv=floor(min([min(kLcfp),min(kLyfp),min(kLrfp)]));
        maxv=ceil(max([max(kLcfp),max(kLyfp),max(kLrfp)]));
        binsize = 0.05;

        bins = minv:binsize:maxv;
        nbins = numel(bins)-1;
        xbin = bins(2:end)-binsize/2;

        [MeanMy, MeanMr, MeanMc, NumMy, NumMr, NumMc, StdEMy, StdEMr, StdEMc, MeanFy, MeanFr, MeanFc, StdEFy, StdEFr, StdEFc] = deal(zeros(nbins,1));
        
        kdifcal = mean(kIcfp);
        
        for j = 1:nbins
            yindex = find(kLyfp>bins(j)&kLyfp<bins(j+1));
            rindex = find(kLrfp>bins(j)&kLrfp<bins(j+1));
            cindex = find(kLcfp>bins(j)&kLcfp<bins(j+1));

            MeanMy(j) = mean(kIyfp(yindex)/ktuscal);
            MeanMr(j) = mean(kIrfp(rindex)/kdnancal);
            MeanMc(j) = mean(kIcfp(cindex)/kdifcal);

            StdEMy(j) = std(kIyfp(yindex)/ktuscal);   %/sqrt(numel(yindex));
            StdEMr(j) = std(kIrfp(rindex)/kdnancal);  %/sqrt(numel(rindex));
            StdEMc(j) = std(kIcfp(cindex)/kdifcal);   %/sqrt(numel(cindex));

            MeanFy(j) = mean(kFyfp(yindex)/ktuscal);
            MeanFr(j) = mean(kFrfp(rindex)/kdnancal);
            MeanFc(j) = mean(kFcfp(cindex)/kdifcal);

            StdEFy(j) = std(kFyfp(yindex)/ktuscal);   %/sqrt(numel(yindex));
            StdEFr(j) = std(kFrfp(rindex)/kdnancal);  %/sqrt(numel(rindex));
            StdEFc(j) = std(kFcfp(cindex)/kdifcal); 
            
            NumMy(j) = numel(yindex);
            NumMr(j) = numel(rindex);
            NumMc(j) = numel(cindex);
        end
        nany = isnan(MeanMy);
        MeanMy(nany)=[]; NumMy(nany)=[]; StdEMy(nany)=[]; MeanFy(nany)=[]; StdEFy(nany)=[]; ybin = xbin(~nany);
        nanr = isnan(MeanMr);
        MeanMr(nanr)=[]; NumMr(nanr)=[]; StdEMr(nanr)=[]; MeanFr(nanr)=[]; StdEFr(nanr)=[]; rbin = xbin(~nanr);
        nanc = isnan(MeanMc);
        MeanMc(nanc)=[]; NumMc(nanc)=[]; StdEMc(nanc)=[]; MeanFc(nanc)=[]; StdEFc(nanc)=[]; cbin = xbin(~nanc);

        % Tus
        figure(33)
        subplot(2,2,2)
        
        hold on
        scatter(kLyfp,kIyfp/ktuscal,30,green,'x','MarkerEdgeAlpha',0.15,'MarkerFaceAlpha',0.15);
        errorbar(ybin,MeanMy,StdEMy,'vertical','LineStyle','-','LineWidth',LW,'Color',green)
%         fity = fit(ybin(1:end-3)',MeanMy(1:end-3),'poly1');
%         fitplot = plot(fity,'g');

        scc = 1;
        h2=scatter(ybin,MeanMy,NumMy*scc,green,'LineWidth',LW);
        scatter(ybin,MeanMy,NumMy*scc,green,'filled','LineWidth',LW,'MarkerFaceAlpha',0.2)
        
        vlplot = vline(bins,'k:');
        set(vlplot,'LineWidth',0.5)

        hold off
        legend off
        xlabel(CC)
        ylabel('Number of proteins')
%         set(fitplot,'LineWidth',2)
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS2,'Color',[0.95 0.95 0.95])
        set(gcf,'Position',pos)
        legend(h2,'Tus')
%         TitlePos = GetTitPos(Tx,Ty);t=title('B.','Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS2)
        
        % DnaN
        figure(33)
        subplot(2,2,1)
        
        hold on
        
        scatter(kLrfp,kIrfp/kdnancal,30,'rx','MarkerEdgeAlpha',0.15,'MarkerFaceAlpha',0.15);
%         scatter(kLrfp,kFrfp/kdnancal,30,oran,'x','MarkerEdgeAlpha',0.15,'MarkerFaceAlpha',0.15);
        errorbar(rbin,MeanMr,StdEMr,'vertical','LineStyle','-','LineWidth',LW,'Color','r')
%         errorbar(rbin,MeanFr,StdEFr,'vertical','LineStyle','-','LineWidth',LW,'Color',oran)
%         fitr = fit(rbin(1:end-3)',MeanMr(1:end-3),'poly1');
%         fitplot = plot(fitr,'r');

        scc = 1.5;
        h3=scatter(rbin,MeanMr,NumMr*scc,'r','LineWidth',LW);
        scatter(rbin,MeanMr,NumMr*scc,'r','filled','LineWidth',LW,'MarkerFaceAlpha',0.2)
        
        vlplot = vline(bins,'k:');
        set(vlplot,'LineWidth',0.5)

        legend off
        hold off
        xlabel(CC)
        ylabel('Number of proteins')
        set(fitplot,'LineWidth',2)
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS2,'Color',[0.95 0.95 0.95])
        set(gcf,'Position',pos)
        legend(h3,'DnaN')
%         TitlePos = GetTitPos(Tx,Ty);t=title('A.','Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS2)
    end
    
  % Numprots2
    if k_numprots2==1
        
        minv=floor(min([min(kLcfp),min(kLyfp),min(kLrfp)]));
        maxv=ceil(max([max(kLcfp),max(kLyfp),max(kLrfp)]));
        binsize = 0.05;

        bins = minv:binsize:maxv;
        nbins = numel(bins)-1;
        xbin = bins(2:end)-binsize/2;

        [MeanMy, MeanMr, MeanMc, NumMy, NumMr, NumMc, StdEMy, StdEMr, StdEMc] = deal(zeros(nbins,1));
        
        kdifcal = mean(kIcfp);
        
        for j = 1:nbins
            yindex = find(kLyfp>bins(j)&kLyfp<bins(j+1));
            rindex = find(kLrfp>bins(j)&kLrfp<bins(j+1));
            cindex = find(kLcfp>bins(j)&kLcfp<bins(j+1));
            
            yind = find(kCLy>bins(j)&kCLy<bins(j+1));
            rind = find(kCLr>bins(j)&kCLr<bins(j+1));
            cind = find(kCLc>bins(j)&kCLc<bins(j+1));
            
            % prot per spot
            MeanMy(j) = mean(kIyfp(yindex)/ktuscal);
            MeanMr(j) = mean(kIrfp(rindex)/kdnancal);
            MeanMc(j) = mean(kIcfp(cindex)/kdifcal);

            StdEMy(j) = std(kIyfp(yindex)/ktuscal); 
            StdEMr(j) = std(kIrfp(rindex)/kdnancal); 
            StdEMc(j) = std(kIcfp(cindex)/kdifcal);

            NumMy(j) = numel(yindex);
            NumMr(j) = numel(rindex);
            NumMc(j) = numel(cindex);
        end
        nany = isnan(MeanMy); MeanMy(nany)=[]; NumMy(nany)=[]; StdEMy(nany)=[]; ybin = xbin(~nany);
        nanr = isnan(MeanMr); MeanMr(nanr)=[]; NumMr(nanr)=[]; StdEMr(nanr)=[]; rbin = xbin(~nanr);
        nanc = isnan(MeanMc); MeanMc(nanc)=[]; NumMc(nanc)=[]; StdEMc(nanc)=[]; cbin = xbin(~nanc);
        
        
        
%         minv2=floor(min([min(kCLc),min(kCLy),min(kCLr)]));
%         maxv2=ceil(max([max(kCLc),max(kCLy),max(kCLr)]));
%         binsize2 = 0.03;
        
        for j = 1:nbins
            yind = find(kCLy>bins(j)&kCLy<bins(j+1));
            rind = find(kCLr>bins(j)&kCLr<bins(j+1));
            cind = find(kCLc>bins(j)&kCLc<bins(j+1));
            
            % Prot per cell
            MIy(j) = mean(kIy(yind)/ktuscal);
            MIr(j) = mean(kIr(rind)/kdnancal);
            MIc(j) = mean(kIc(cind)/kdifcal);

            StdIy(j) = std(kIy(yind)/ktuscal); 
            StdIr(j) = std(kIr(rind)/kdnancal); 
            StdIc(j) = std(kIc(cind)/kdifcal);
            
            % Full cell
            MFy(j) = mean(kFy(yind)/ktuscal);
            MFr(j) = mean(kFr(rind)/kdnancal);
            MFc(j) = mean(kFc(cind)/kdifcal);

            StdFy(j) = std(kFy(yind)/ktuscal); 
            StdFr(j) = std(kFr(rind)/kdnancal); 
            StdFc(j) = std(kFc(cind)/kdifcal);
            
            % Quotient
            MQy(j) = mean(kQy(yind));
            MQr(j) = mean(kQr(rind));
            MQc(j) = mean(kQc(cind));

            StdQy(j) = std(kQy(yind)); 
            StdQr(j) = std(kQr(rind)); 
            StdQc(j) = std(kQc(cind));
            
            NMy(j) = numel(yind);
            NMr(j) = numel(rind);
            NMc(j) = numel(cind);
        end
        ny = isnan(MIy); MIy(ny)=[]; NMy(ny)=[]; StdIy(ny)=[]; MFy(ny)=[]; StdFy(ny)=[]; MQy(ny)=[]; StdQy(ny)=[]; ybi = xbin(~ny);
        nr = isnan(MIr); MIr(nr)=[]; NMr(nr)=[]; StdIr(nr)=[]; MFr(nr)=[]; StdFr(nr)=[]; MQr(nr)=[]; StdQr(nr)=[]; rbi = xbin(~nr);
        nc = isnan(MIc); MIc(nc)=[]; NMc(nc)=[]; StdIc(nc)=[]; MFc(nc)=[]; StdFc(nc)=[]; MQc(nc)=[]; StdQc(nc)=[]; cbi = xbin(~nc);
        
        % DnaN Per Spot
        figure(44)
        set(gcf,'Position',screensize3)
        subplot(2,2,1)
        
        hold on
%         scatter(kLrfp,kIrfp/kdnancal,30,'rx','MarkerEdgeAlpha',0.15,'MarkerFaceAlpha',0.15);
%         errorbar(rbin,MeanMr,StdEMr,'vertical','LineStyle','-','LineWidth',LW,'Color','r')
% 
%         scc = 1.5;
%         h3=scatter(rbin,MeanMr,NumMr*scc,'r','LineWidth',LW);
%         scatter(rbin,MeanMr,NumMr*scc,'r','filled','LineWidth',LW,'MarkerFaceAlpha',0.2)

        X1=[rbin,fliplr(rbin)];
        Y1=[(MeanMr-StdEMr)',fliplr((MeanMr+StdEMr)')];
        fi = fill(X1,Y1,red);
        set(fi,'LineStyle','none','FaceAlpha',0.2)
        plot(rbin,(MeanMr+StdEMr)','-','LineWidth',LW,'Color',red);
        plot(rbin,(MeanMr-StdEMr)','-','LineWidth',LW,'Color',red);
        plot(rbin,MeanMr','-','LineWidth',3*LW,'Color',red);
        
        axis([0 1 0 5])
        vl = vline(0.9,'k--');
        set(vl,'LineWidth',LW)
        hold off
        legend off
        xlabel(CC)
        ylabel('Proteins per spot')
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS2,'YGrid','on','GridLineStyle','-','Color',[0.95 0.95 0.95])
        TitlePos = GetTitPos(Tx,Ty);t=title('A.','Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS2)
        
        % DnaN Per Cell
        figure(44)
        subplot(2,2,2)
        
        hold on
        X2=[rbi,fliplr(rbi)];
        Y1=[MIr-StdIr,fliplr(MIr+StdIr)];
        fila = fill(X2,Y1,'r');
        set(fila,'LineStyle','none','FaceAlpha','0.1')
        plot(rbi,MIr+StdIr,'r-','LineWidth',LW);
        plot(rbi,MIr-StdIr,'r-','LineWidth',LW);
        plot(rbi,MIr,'-','LineWidth',3*LW,'Color','r')
        
        axis([0 1 -2 7])
        vl = vline(0.9,'k--');
        set(vl,'LineWidth',LW)
        hold off
        legend off
        xlabel(CC)
        ylabel('DNA-bound proteins per cell')
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS2,'YGrid','on','GridLineStyle','-','Color',[0.95 0.95 0.95])
        TitlePos = GetTitPos(Tx,Ty);t=title('B.','Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS2)

        % DnaN Per Cell % FCI
        figure(44)
        subplot(2,2,3)
        
        hold on
%         scatter(kLrfp,kIrfp/kdnancal,30,'rx','MarkerEdgeAlpha',0.15,'MarkerFaceAlpha',0.15);
        X2=[rbi,fliplr(rbi)];
        Y1=[MIr-StdIr,fliplr(MIr+StdIr)];
        fila = fill(X2,Y1,'r');
        set(fila,'LineStyle','none','FaceAlpha','0.1')
        plot(rbi,MIr+StdIr,'r-','LineWidth',LW);
        plot(rbi,MIr-StdIr,'r-','LineWidth',LW);
        h1=plot(rbi,MIr,'-','LineWidth',3*LW,'Color','r');
        
%         scatter(kLrfp,kFrfp/kdnancal,30,blue,'x','MarkerEdgeAlpha',0.15,'MarkerFaceAlpha',0.15);
        Y2=[MFr-StdFr,fliplr(MFr+StdFr)];
        filb = fill(X2,Y2,blue);
        set(filb,'LineStyle','none','FaceAlpha','0.1')
        plot(rbi,MFr+StdFr,'-','LineWidth',LW,'Color',blue);
        plot(rbi,MFr-StdFr,'-','LineWidth',LW,'Color',blue);
        h2=plot(rbi,MFr,'-','LineWidth',3*LW,'Color',blue);
        vl = vline(0.9,'k--');
        set(vl,'LineWidth',LW)
        hold off
        legend off
        xlabel(CC)
        ylabel('Proteins per cell')
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS2,'YGrid','on','GridLineStyle','-','Color',[0.95 0.95 0.95])
        legend([h1,h2],{'DNA-bound','Total'},'Location','Northwest')
        TitlePos = GetTitPos(Tx,Ty);t=title('C.','Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS2)
        clear t
        
        % Persentage
        figure(44)
        subplot(2,2,4)
        
        hold on
        X2=[rbi,fliplr(rbi)];
        Y1=[MQr-StdQr,fliplr(MQr+StdQr)];

        fila = fill(X2,Y1,oran);
        set(fila,'LineStyle','none','FaceAlpha','0.1')
        plot(rbi,MQr-StdQr,'-','LineWidth',LW,'Color',oran);
        plot(rbi,MQr+StdQr,'-','LineWidth',LW,'Color',oran);
        plot(rbi,MQr,'-','LineWidth',3*LW,'Color',oran)
        
        axis([0,1,0,25])
        vl = vline(0.9,'k--');
        set(vl,'LineWidth',LW)
        hold off
        legend off
        xlabel(CC)
        ylabel('DNA-bound DnaN (%)')
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS2,'YGrid','on','GridLineStyle','-','Color',[0.95 0.95 0.95])
        TitlePos = GetTitPos(Tx,Ty);t=title('D.','Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS2)

        
        
        % Tus
        % Tus Per Spot
        figure(55)
        subplot(2,2,1)
        
        hold on
%         scatter(kLrfp,kIrfp/kdnancal,30,'rx','MarkerEdgeAlpha',0.15,'MarkerFaceAlpha',0.15);
%         errorbar(ybin,MeanMy,StdEMy,'vertical','LineStyle','-','LineWidth',LW,'Color','r')
% 
%         scc = 1.5;
%         h3=scatter(ybin,MeanMy,NumMr*scc,'r','LineWidth',LW);
%         scatter(ybin,MeanMy,NumMr*scc,'r','filled','LineWidth',LW,'MarkerFaceAlpha',0.2)

        X1=[ybin,fliplr(ybin)];
        Y1=[(MeanMy-StdEMy)',fliplr((MeanMy+StdEMy)')];
        fi = fill(X1,Y1,[0,0.2,0.05]);
        set(fi,'LineStyle','none','FaceAlpha',0.2)
        plot(ybin,(MeanMy+StdEMy)','-','LineWidth',LW,'Color',[0,0.2,0.05]);
        plot(ybin,(MeanMy-StdEMy)','-','LineWidth',LW,'Color',[0,0.2,0.05]);
        plot(ybin,MeanMy','-','LineWidth',3*LW,'Color',[0,0.2,0.05]);
        
        axis([0 1 0 5])
        vl = vline(0.9,'k--');
        set(vl,'LineWidth',LW)
        hold off
        legend off
        xlabel(CC)
        ylabel('Proteins per spot')
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS2,'YGrid','on','GridLineStyle','-','Color',[0.95 0.95 0.95])
        TitlePos = GetTitPos(Tx,Ty);t=title('A.','Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS2)
        
        % Tus Per Cell
        figure(55)
        subplot(2,2,2)
        
        hold on
        X2=[ybi,fliplr(ybi)];
        Y1=[MIy-StdIy,fliplr(MIy+StdIy)];
        fila = fill(X2,Y1,green);
        set(fila,'LineStyle','none','FaceAlpha','0.1')
        plot(ybi,MIy+StdIy,'-','LineWidth',LW,'Color',green);
        plot(ybi,MIy-StdIy,'-','LineWidth',LW,'Color',green);
        plot(ybi,MIy,'-','LineWidth',3*LW,'Color',green')
        
        axis([0 1 -2 7])
        vl = vline(0.9,'k--');
        set(vl,'LineWidth',LW)
        hold off
        legend off
        xlabel(CC)
        ylabel('DNA-bound proteins per cell')
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS2,'YGrid','on','GridLineStyle','-','Color',[0.95 0.95 0.95])
        TitlePos = GetTitPos(Tx,Ty);t=title('B.','Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS2)

        % Tus Per Cell % FCI
        figure(55)
        subplot(2,2,3)
        
        hold on
%         scatter(kLrfp,kIrfp/kdnancal,30,'rx','MarkerEdgeAlpha',0.15,'MarkerFaceAlpha',0.15);
        X2=[ybi,fliplr(ybi)];
        Y1=[MIy-StdIy,fliplr(MIy+StdIy)];
        fila = fill(X2,Y1,green);
        set(fila,'LineStyle','none','FaceAlpha','0.1')
        plot(ybi,MIy+StdIy,'-','LineWidth',LW,'Color',green);
        plot(ybi,MIy-StdIy,'-','LineWidth',LW,'Color',green);
        h1=plot(ybi,MIy,'-','LineWidth',3*LW,'Color',green);
        
%         scatter(kLrfp,kFrfp/kdnancal,30,blue,'x','MarkerEdgeAlpha',0.15,'MarkerFaceAlpha',0.15);
        Y2=[MFy-StdFy,fliplr(MFy+StdFy)];
        filb = fill(X2,Y2,blue);
        set(filb,'LineStyle','none','FaceAlpha','0.1')
        plot(ybi,MFy+StdFy,'-','LineWidth',LW,'Color',blue);
        plot(ybi,MFy-StdFy,'-','LineWidth',LW,'Color',blue);
        h2=plot(ybi,MFy,'-','LineWidth',3*LW,'Color',blue);
        vl = vline(0.9,'k--');
        set(vl,'LineWidth',LW)
        hold off
        legend off
        xlabel(CC)
        ylabel('Proteins per cell')
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS2,'YGrid','on','GridLineStyle','-','Color',[0.95 0.95 0.95])
        set(gcf,'Position',pos)
        legend([h1,h2],{'DNA-bound','Total'},'Location','Northwest')
        TitlePos = GetTitPos(Tx,Ty);t=title('C.','Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS2)
        
        % Persentage
        figure(55)
        subplot(2,2,4)
        
        hold on
        X2=[ybi,fliplr(ybi)];
        Y1=[MQy-StdQy,fliplr(MQy+StdQy)];

        fila = fill(X2,Y1,oran);
        set(fila,'LineStyle','none','FaceAlpha','0.1')
        plot(ybi,MQy-StdQy,'-','LineWidth',LW,'Color',oran);
        plot(ybi,MQy+StdQy,'-','LineWidth',LW,'Color',oran);
        plot(ybi,MQy,'-','LineWidth',3*LW,'Color',oran)
        
        axis([0,1,0,25])
        vl = vline(0.9,'k--');
        set(vl,'LineWidth',LW)
        hold off
        legend off
        xlabel(CC)
        ylabel('DNA-bound Tus (%)')
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS2,'YGrid','on','GridLineStyle','-','Color',[0.95 0.95 0.95])
        TitlePos = GetTitPos(Tx,Ty);t=title('D.','Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS2)

        set(gcf,'Position',screensize3)
    end
    
%% Agar
    % Get Data
    if a_getdata==1    
        a_channels=[1 2 3 4 5 7 8 9];
        aIntensityval = [700, 300, 700];
        umperpx=0.159;

        aLcfp=[];    aLyfp=[];    aLrfp=[];
        aPcfp=[];    aPyfp=[];    aPrfp=[];
        aIcfp=[];    aIyfp=[];    aIrfp=[];
        aFcfp=[];    aFyfp=[];    aFrfp=[];
        aNSc=[];     aNSy=[];     aNSr=[];
        aNSc=[];     aNSy=[];     aNSr=[];
        aCL = [];

        yfp.filterval=8000;

        j=1; 
        for i=a_channels
            E{j}=load(strcat(folder,slash,'Results',num2str(i),'.mat')); 
            j=j+1;
        end

        allCFP_L = [];

        Nexp=size(E,2);
        
        % i: exp number
        % j: cell number
        % k: spot number

        for i=1:Nexp

            Ncells{i}=size(E{i}.DataStruct,2);


            for j=1:Ncells{i} 

                if ~isempty(E{i}.DataStruct(1,j).Lnorm)        
                    LNormCFP{i,j}=E{i}.DataStruct(1,j).Lnorm;
                else
                    LNormCFP{i,j}=0;
                end

                CellLength{i,j}=E{i}.DataStruct(1,j).CellLength;
                aCL = [aCL CellLength{i,j}];

                CFPld{i,j}=E{i}.DataStruct(1,j).ld;
                YFPld{i,j}=E{i}.DataStruct(2,j).ld;
                RFPld{i,j}=E{i}.DataStruct(3,j).ld;

                NspotsCFP=size(CFPld{i,j},2);
                NspotsYFP=size(YFPld{i,j},2);
                NspotsRFP=size(RFPld{i,j},2);        

                if NspotsCFP==0
                    CFPld{i,j}{1}=[];
                else
                    for k=1:NspotsCFP
                        aLcfp=[aLcfp CellLength{i,j}];
                        aPcfp=[aPcfp CFPld{i,j}{k}(1,2)/CellLength{i,j}];
                        aIcfp=[aIcfp CFPld{i,j}{k}(1,1)];
                        aFcfp=[aFcfp CFPld{i,j}{k}(1,7)];
                    end
                end
                aNSc = [aNSc NspotsCFP];

                if NspotsYFP==0
                    YFPld{i,j}{1}=[];
                else
                    for k=1:NspotsYFP
                        aLyfp=[aLyfp CellLength{i,j}];
                        aPyfp=[aPyfp YFPld{i,j}{k}(1,2)/CellLength{i,j}];
                        aIyfp=[aIyfp YFPld{i,j}{k}(1,1)];
                        aFyfp=[aFyfp YFPld{i,j}{k}(1,7)];
                    end
                end
                aNSy = [aNSy NspotsYFP];

                if NspotsRFP==0
                    RFPld{i,j}{1}=[];
                else
                    for k=1:NspotsRFP            
                        aLrfp=[aLrfp CellLength{i,j}];
                        aPrfp=[aPrfp RFPld{i,j}{k}(1,2)/CellLength{i,j}];
                        aIrfp=[aIrfp RFPld{i,j}{k}(1,1)];
                        aFrfp=[aFrfp RFPld{i,j}{k}(1,7)];
                    end
                end
                aNSr = [aNSr NspotsRFP];
            end
        end

        % tozero = Iyfp<yfp.filterval;
        % Iyfp(tozero) = 0;
        % Lyfp(tozero) = 0;
        % Pyfp(tozero) = 0;
        % Fyfp(tozero) = 0;
    end

    % Position vs. cell length
    if a_PvsCL==1

        figure
        set(gcf, 'position', screensize)

        % CFP
        subplot(1,3,1)
        hold on
        scatter(single(aLcfp),aPcfp,aIcfp,blue,'filled');
        myfit=polyfit(aLcfp,aPcfp,4);
        x=15:0.1:45;
        y=polyval(myfit,x);
        %plot(x,y,'r','LineWidth',5)
        ylabel(P); 
        t=title('Agar data: CFP');
        hold off
        axis([CLbound 0 1])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font)


        % YFP
        subplot(1,3,2);
        hold on
        scatter(single(aLyfp),aPyfp,aIyfp,green,'filled');
        myfit=polyfit(aLyfp,aPyfp,4);
        x=15:0.1:45;
        y=polyval(myfit,x);
        % plot(x,y,'k','LineWidth',5)
        xlabel(CL);
        title('Agar data: YFP')
        hold off
        axis([CLbound 0 1])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font)

        % RFP
        subplot(1,3,3)
        hold on
        scatter(single(aLrfp),aPrfp,aIrfp,'r','filled');
        myfit=polyfit(aLrfp,aPrfp,4);
        x=15:0.1:45;
        y=polyval(myfit,x);
        title('Agar data: RFP')
        % plot(x,y,'k','LineWidth',5)
        hold off
        axis([CLbound 0 1])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font)
    end

    % Intensity vs. position
    if a_IvsP==1
        
        figure
        set(gcf, 'position', screensize)
        
        % CFP
        subplot(1,3,1)
        hold on
        scatter(aPcfp,aIcfp,'x','MarkerEdgeColor',blue);
        myfit=polyfit(aPcfp,aIcfp,4);
        x=0:0.001:1;
        y=polyval(myfit,x);
        plot(x,y,'k','LineWidth',3)
        ylabel(SI); 
        t=title('Agar data: CFP');
        hold off
        axis([0 1 -0.1 90])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font)


        % YFP
        subplot(1,3,2)
        hold on
        scatter(aPyfp,aIyfp,'x','MarkerEdgeColor',green);
        myfit=polyfit(aPyfp,aIyfp,4);
        x=0:00.1:1;
        y=polyval(myfit,x);
        plot(x,y,'k','LineWidth',3)
        xlabel(P);
        t=title('Agar data: YFP');
        hold off
        axis([0 1 -0.1 90])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font)


        % RFP
        subplot(1,3,3)
        hold on
        scatter(aPrfp,aIrfp,'r','x');
        myfit=polyfit(aPrfp,aIrfp,4);
        x=0:00.1:1;
        y=polyval(myfit,x);
        plot(x,y,'k','LineWidth',3)
        t=title('Agar data: RFP');
        hold off
        axis([0 1 -0.1 90])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font)
    end

    % Numspots vs. position
    if a_NSvsP==1
        bins = 15;
        thisedge = (0:bins)/bins;

        figure
        set(gcf, 'position', screensize)

        % CFP
        subplot(1,3,1)
        [numbin,edges] = histcounts(aPcfp,thisedge);
        norm = max(numbin)/80;
        X = diff(edges);
        X = cumsum(X) - X(1)/2;
        hold on
        scatter(aPcfp,aIcfp,'x','MarkerEdgeColor',blue);
        plot(X,numbin/norm,'k','LineWidth',3)
        hold off
        ylabel(SIN);
        t=title('Agar data: CFP');
        axis([0 1 -0.1 90])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font)

        % YFP
        subplot(1,3,2)
        [numbin,edges] = histcounts(aPyfp,thisedge);
        norm = max(numbin)/80;
        X = diff(edges);
        X = cumsum(X) - X(1)/2;
        hold on
        scatter(aPyfp,aIyfp,'x','MarkerEdgeColor',green);
        plot(X,numbin/norm,'k','LineWidth',3)
        hold off
        xlabel(P);
        t=title('Agar data: YFP');
        axis([0 1 -0.1 90])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font)

        % RFP
        subplot(1,3,3)
        [numbin,edges] = histcounts(aPrfp,thisedge);
        norm = max(numbin)/80;
        X = diff(edges);
        X = cumsum(X) - X(1)/2;
        hold on
        scatter(aPrfp,aIrfp,'r','x');
        plot(X,numbin/norm,'k','LineWidth',3)
        hold off
        t=title('Agar data: RFP');
        axis([0 1 -0.1 90])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font)
    end

    % Numspots vs. position & cell length
    if a_Heatmap==1
        
        bins = 15;
        thisedge2{1} = linspace(min(aLcfp),max(aLcfp),bins+5);
        thisedge2{2} = (0:bins)/bins;

        % CFP

        subplot(1,3,1)
        aNumcfp(1,:) = aLcfp;
        aNumcfp(2,:) = aPcfp;

        %Filter on Intensity
        [bin_cfp,idx_cfp]=find(kIcfp*aIntensityval(1)>6000); 

        NumCFP(:,idx_cfp)=aNumcfp(:,idx_cfp);
        NumCFPnz(1,:)=nonzeros(NumCFP(1,:));
        NumCFPnz(2,:)=nonzeros(NumCFP(2,:));

        FilteredSignalCFP=size(NumCFPnz,2)/size(kIcfp,2); % Percentage of total signal

        figure(21)
        subplot(3,2,2)
        Heatmap = hist3(NumCFPnz','Edges',thisedge2);
        pcolor(thisedge2{1},(thisedge2{2}),Heatmap');
        colormap(hot) % heat map
        t=title('Agar data: CFP');
        grid on
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font)

%         figure(32)
%         subplot(1,3,1)
%         hold on
%         hist3(NumCFPnz','Edges',thisedge2)
%         colormap(winter) % heat map
%         set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
%         xlabel(CL); ylabel(P); zlabel(NS)
%         t=title('Agar data: CFP');
%         grid off
%         hold off
%         axis([min(aLcfp),max(aLcfp),0,1])
%         view(3)
%         set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
%         set(t,'Fontname',Font)

        % YFP
        aNumyfp(1,:) = aLyfp;
        aNumyfp(2,:) = aPyfp;

        %Filter on Intensity
        [bin_yfp,idx_yfp]=find(kIyfp*aIntensityval(2)>yfp.filterval); 

        aNumYFP(:,idx_yfp)=aNumyfp(:,idx_yfp);
        aNumYFPnz(1,:)=nonzeros(aNumYFP(1,:));
        aNumYFPnz(2,:)=nonzeros(aNumYFP(2,:));

        FilteredSignalYFP=size(aNumYFPnz,2)/size(kIyfp,2); % Percentage of total signal

        figure(21)
        subplot(3,2,4)
        Heatmap = hist3(aNumYFPnz','Edges',thisedge2);
        pcolor(thisedge2{1},(thisedge2{2}),Heatmap');
        colormap(hot) % heat map
        t=title('Agar data: YFP');
        grid on
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font)

%         figure(32)
%         subplot(1,3,2)
%         hold on
%         hist3(aNumYFPnz','Edges',thisedge2)
%         colormap(autumn) % heat map
%         set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
%         xlabel(CC); ylabel(P); 
%         t=title('Agar data: YFP');
%         grid off
%         hold off
%         axis([min(aLyfp),max(aLyfp),0,1])
%         view(3)
%         set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
%         set(t,'Fontname',Font)

        % RFP
        aNumrfp(1,:) = aLrfp;
        aNumrfp(2,:) = aPrfp;

        %Filter on Intensity
        [bin_rfp,idx_rfp]=find(kIrfp*aIntensityval(3)>10000); 

        NumRFP(:,idx_rfp)=aNumrfp(:,idx_rfp);

        NumRFPnz(1,:)=nonzeros(NumRFP(1,:));
        NumRFPnz(2,:)=nonzeros(NumRFP(2,:));

        FilteredSignalRFP=size(NumRFPnz,2)/size(kIrfp,2);

        figure(21)
        subplot(3,2,6)
        Heatmap = hist3(NumRFPnz','Edges',thisedge2);
        h = pcolor(thisedge2{1},(thisedge2{2}),Heatmap');
        colormap(hot) % heat map
        xlabel(CL);
        t=title('Agar data: RFP');
        grid on
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(t,'Fontname',Font)


%         figure(32)
%         subplot(1,3,3)
%         hold on
%         hist3(NumRFPnz','Edges',thisedge2)
%         colormap(hot) % heat map
%         set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
%         xlabel(CL); ylabel(P); 
%         t=title('Agar data: RFP');
%         grid off
%         hold off
%         axis([min(aLrfp),max(aLrfp),0,1])
%         view(3)
%         set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
%         set(t,'Fontname',Font)

    end

    % Full cell intensity vs. celllength
    if a_FCIvsCL==1
        clear plotcfp plotyfp plotrfp
        fig6 = figure(6);
        set(fig6,'Position',[20,300,1800,500])

        % CFP

        plotcfp(1,:) = aLcfp;
        plotcfp(2,:) = aFcfp;
        plotcfp(3,:) = aIcfp*aIntensityval(1);

        plotcfp = unique(plotcfp','rows')';

        subplot(1,3,1)
        hold on
        scatter(plotcfp(1,:),plotcfp(2,:),'b','o','filled');
        scatter(plotcfp(1,:),plotcfp(3,:),'r','o','filled');
        myfit=polyfit(plotcfp(1,:),plotcfp(2,:),1);
        myfit2=polyfit(plotcfp(1,:),plotcfp(3,:),1);
        x=12:0.1:43;
        y=polyval(myfit,x);
        y2=polyval(myfit2,x);
        plot(x,y,'b','LineWidth',3)
        plot(x,y2,'r','LineWidth',3)
        xlabel('Cell Length (px)'); ylabel('Intensity (-)'); 
        title('Agar data: CFP')
        hold off
        axis([12 43 -0.1 2*10^5])
        set(gca,'FontSize',16)
        legend('Full cell intensity','Spot intensity')


        % YFP

        plotyfp(1,:) = aLyfp;
        plotyfp(2,:) = aFyfp;
        plotyfp(3,:) = aIyfp*aIntensityval(2);

        plotyfp = unique(plotyfp','rows')';

        figure(1)
        hold on
        scatter(plotyfp(1,:),plotyfp(2,:),'b','o','filled');
        scatter(plotyfp(1,:),plotyfp(3,:),'r','o','filled');
        myfit=polyfit(plotyfp(1,:),plotyfp(2,:),1);
        myfit2=polyfit(plotyfp(1,:),plotyfp(3,:),1);
        x=12:0.1:43;
        y=polyval(myfit,x);
        y2=polyval(myfit2,x);
        plot(x,y,'b','LineWidth',3)
        % plot(x,y2,'r','LineWidth',3)
        xlabel('Cell Length'); ylabel('Intensity'); 
        title('Tus Stoichiometry vs. Length')
        hold off
        axis([12 43 -0.1 2*10^5])
        set(gca,'FontSize',16)


        % RFP

        plotrfp(1,:) = aLrfp;
        plotrfp(2,:) = aFrfp;
        plotrfp(3,:) = aIrfp*aIntensityval(3);
        plotrfp = unique(plotrfp','rows')';

        subplot(1,3,3)
        hold on
        scatter(plotrfp(1,:),plotrfp(2,:),'b','o','filled');
        scatter(plotrfp(1,:),plotrfp(3,:),'r','o','filled');
        myfit=polyfit(plotrfp(1,:),plotrfp(2,:),1);
        myfit2=polyfit(plotrfp(1,:),plotrfp(3,:),1);
        x=12:0.1:43;
        y=polyval(myfit,x);
        y2=polyval(myfit2,x);
        plot(x,y,'b','LineWidth',3)
        plot(x,y2,'r','LineWidth',3)
        xlabel('Cell Length'); ylabel('Normalized full cell intensity'); 
        title('Agar data: RFP')
        hold off
        axis([12 43 -0.1 4*10^5])
        set(gca,'FontSize',16)
    end
        
    % Numspots/cell vs. cell length
    if a_royplots==1
        fig1 = figure(7);

        for i=1:max(aNSr);
        p{i} = find(aNSr == i);
        m(i) = mean(aCL(p{i}));
        Spotstd(i)=std(aCL(p{i}));
        end

        myfit=polyfit(m,1:4,1);
        x=[12,43];
        y=polyval(myfit,x);

        Y=linspace(1,max(aNSr),max(aNSr));

        Xl=m-Spotstd;
        Xr=m+Spotstd;

        hold on
        scatter(aCL,aNSr,'m')
        plot(m,Y,'b--','LineWidth',3)
        plot(Xl,Y,'b','LineWidth',1)
        plot(Xr,Y,'b','LineWidth',1)
        axis([12, 43, 0, 5])
        xlabel('Cell length'); ylabel('Number of spots per cell')
        title('Numspots/cell vs. cell length for YFP')
        hold off
        set(gca,'FontSize',16)
        legend('Spot vs Length points','Mean value','std')
    end 
    
    % Numspots
    if a_numspots ==1
        minv=floor(min(aCL));
        maxv=ceil(max(aCL));
        binsize = 2;

        bins = minv:binsize:maxv;
        nbins = numel(bins)-1;
        xbin = bins(2:end)-binsize/2;

        [MeanMy, MeanMr, MeanMc, NumMy, NumMr, NumMc, StdEMy, StdEMr, StdEMc] = deal(zeros(nbins,1));

        for j = 1:nbins
            aindex = find(aCL>bins(j)&aCL<bins(j+1));

            thismy = mean(aNSy(aindex));
            thismr = mean(aNSr(aindex));
            thismc = mean(aNSc(aindex));

            MeanMy(j) = thismy;
            MeanMr(j) = thismr;
            MeanMc(j) = thismc;

            StdEMy(j) = std(aNSy(aindex)); %/sqrt(numel(aindex));
            StdEMr(j) = std(aNSr(aindex)); %/sqrt(numel(aindex));
            StdEMc(j) = std(aNSc(aindex)); %/sqrt(numel(aindex));

            NumMy(j) = numel(aindex);
            NumMr(j) = numel(aindex);
            NumMc(j) = numel(aindex);
        end
        nany = isnan(MeanMy);
        MeanMy(nany)=[]; NumMy(nany)=[]; StdEMy(nany)=[]; ybin = xbin(~nany);
        nanr = isnan(MeanMr);
        MeanMr(nanr)=[]; NumMr(nanr)=[]; StdEMr(nanr)=[]; rbin = xbin(~nanr);
        nanc = isnan(MeanMc);
        MeanMc(nanc)=[]; NumMc(nanc)=[]; StdEMc(nanc)=[]; cbin = xbin(~nanc);

        fity = fit(ybin(1:end-3)'*px,MeanMy(1:end-3),'poly1');
        fitr = fit(rbin(1:end-3)'*px,MeanMr(1:end-3),'poly1');
        fitc = fit(cbin'*px,MeanMc,'poly1');

        figure
        hold on
        
        shift  = 0.02;
        
        errorbar(rbin*px+shift,MeanMr,StdEMr,'vertical','LineStyle','none','LineWidth',LW,'Color',[1,0,0])
        errorbar(ybin*px,MeanMy,StdEMy,'vertical','LineStyle','none','LineWidth',LW,'Color',green)
        errorbar(cbin*px-shift,MeanMc,StdEMc,'vertical','LineStyle','none','LineWidth',LW,'Color',blue)
        
        b=scatter(aCL*px,aNSr,30,'rx','MarkerEdgeAlpha',0.2,'MarkerFaceAlpha',0.2);
        a=scatter(aCL*px,aNSy,30,green,'x','MarkerEdgeAlpha',0.15,'MarkerFaceAlpha',0.15);
        c=scatter(aCL*px,aNSc,30,blue,'x','MarkerEdgeAlpha',0.1,'MarkerFaceAlpha',0.1);

        % plot(fitc,'b')
        plotfit1 = plot(fity,'g');
        plotfit2 = plot(fitr,'r');
        legend off
        axis([10*px 45*px 0 6])

        scc = 2;
        h1=scatter(cbin*px-shift,MeanMc,NumMc*scc,blue,'LineWidth',LW);
        h2=scatter(ybin*px,MeanMy,NumMy*scc,green,'LineWidth',LW);
        h3=scatter(rbin*px+shift,MeanMr,NumMr*scc,'r','LineWidth',LW);
        scatter(cbin*px-shift,MeanMc,NumMc*scc,blue,'filled','LineWidth',LW,'MarkerFaceAlpha',0.2);
        scatter(ybin*px,MeanMy,NumMy*scc,green,'filled','LineWidth',LW,'MarkerFaceAlpha',0.2);
        scatter(rbin*px+shift,MeanMr,NumMr*scc,'r','filled','LineWidth',LW,'MarkerFaceAlpha',0.2);

        vl1=vline(27*px,'k--');
        vl2=vline(33*px,'k--');
        
        vlplot = vline(bins*px,'k:');
        set(vlplot,'LineWidth',0.5)
        
        hold off

        xlabel('Cell length (\mum)')
        ylabel('Number of spots')

        txt1=text(20,5.5,'I','HorizontalAlignment','center');
        txt2=text(30,5.5,'II','HorizontalAlignment','center');
        txt3=text(40,5.5,'III','HorizontalAlignment','center');
        set(txt1,'Fontname','Raleway Black','FontSize',FS)
        set(txt2,'Fontname','Raleway Black','FontSize',FS)
        set(txt3,'Fontname','Raleway Black','FontSize',FS)

        set(vl1,'LineWidth',2)
        set(vl2,'LineWidth',2)
        set(plotfit1,'LineWidth',2)
        set(plotfit2,'LineWidth',2)
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        set(gcf,'Position',pos)
        legend([h1,h2,h3],{'dif','Tus','DnaN'})
        set(gca,'YGrid','on')
        set(gca,'GridLineStyle','-')
        
        text(20*px,5.5,'I','Fontname',TitleFont,'FontSize',TitleFS)
        text(30*px,5.5,'II','Fontname',TitleFont,'FontSize',TitleFS)
        text(40*px,5.5,'III','Fontname',TitleFont,'FontSize',TitleFS)
        
        if true
            figure
            hold on
            ybin = ybin(1:end-3);
            MeanMy = MeanMy(1:end-3);
            StdEMy = StdEMy(1:end-3);
            
            X2=[ybin*px,fliplr(ybin*px)];
            Y1=[(MeanMy-StdEMy)',fliplr((MeanMy+StdEMy)')];
            fila = fill(X2,Y1,green);
            set(fila,'LineStyle','none','FaceAlpha','0.1')
            plot(ybin*px,MeanMy+StdEMy,'-','LineWidth',LW,'Color',green);
            plot(ybin*px,MeanMy-StdEMy,'-','LineWidth',LW,'Color',green);
            plot(ybin*px,MeanMy,'-','LineWidth',3*LW,'Color',green)        

            axis([10*px 45*px-1 0 6])
            vl=vline(27*px,'k--');
            set(vl,'LineWidth',LW)
            hold off
            legend off
            xlabel(CL)
            ylabel('Proteins per cell')
            set(gcf,'Position',pos)
            set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS,'YGrid','on','GridLineStyle','-')
        end

    end
    
    % Numprots
    if a_numprots==1
        minv=floor(min([min(aLcfp),min(aLyfp),min(aLrfp)]*px));
        maxv=ceil(max([max(aLcfp),max(aLyfp),max(aLrfp)]*px));
        binsize = 0.075*3;
        
        maxv=4.4;
        binsize = 0.1;

        bins = minv:binsize:maxv;
        nbins = numel(bins)-1;
        xbin = bins(2:end)-binsize/2;

        [MeanMy, MeanMr, MeanMc, NumMy, NumMr, NumMc, StdEMy, StdEMr, StdEMc] = deal(zeros(nbins,1));
        
        adifcal = mean(aIcfp);
        
        for j = 1:nbins
            yindex = (aLyfp*px>bins(j)&aLyfp*px<bins(j+1));
            rindex = (aLrfp*px>bins(j)&aLrfp*px<bins(j+1));
            cindex = (aLcfp*px>bins(j)&aLcfp*px<bins(j+1));
            
            NumMy(j) = sum(yindex);
            NumMr(j) = sum(rindex);
            NumMc(j) = sum(cindex);

            MeanMy(j) = mean(aIyfp(yindex)/atuscal);
            MeanMr(j) = mean(aIrfp(rindex)/adnancal);
            MeanMc(j) = mean(aIcfp(cindex)/adifcal);

            StdEMy(j) = std(aIyfp(yindex)/atuscal);     %/sqrt(NumMy(j));
            StdEMr(j) = std(aIrfp(rindex)/adnancal);    %/sqrt(NumMr(j));
            StdEMc(j) = std(aIcfp(cindex)/adifcal);     %/sqrt(NumMc(j));
            
            if true % remove outliers
                TIyfp = aIyfp(yindex)/atuscal;
                TIrfp = aIrfp(rindex)/adnancal;
                TIcfp = aIcfp(cindex)/adifcal;
                
                yi2 = (TIyfp>MeanMy(j)-StdEMy(j)&TIyfp<MeanMy(j)+StdEMy(j));
                ri2 = (TIrfp>MeanMr(j)-StdEMr(j)&TIrfp<MeanMr(j)+StdEMr(j));
                ci2 = (TIcfp>MeanMc(j)-StdEMc(j)&TIcfp<MeanMc(j)+StdEMc(j));
                
                NumMy(j) = sum(yi2);
                NumMr(j) = sum(ri2);
                NumMc(j) = sum(ci2);

                MeanMy(j) = mean(TIyfp(yi2));
                MeanMr(j) = mean(TIrfp(ri2));
                MeanMc(j) = mean(TIcfp(ci2));

                StdEMy(j) = std(TIyfp(yi2));
                StdEMr(j) = std(TIrfp(ri2));
                StdEMc(j) = std(TIcfp(ci2));
            end
            
        end
        nany = isnan(MeanMy);
        MeanMy(nany)=[]; NumMy(nany)=[]; StdEMy(nany)=[]; ybin = xbin(~nany);
        nanr = isnan(MeanMr);
        MeanMr(nanr)=[]; NumMr(nanr)=[]; StdEMr(nanr)=[]; rbin = xbin(~nanr);
        nanc = isnan(MeanMc);
        MeanMc(nanc)=[]; NumMc(nanc)=[]; StdEMc(nanc)=[]; cbin = xbin(~nanc);

        % Tus
        figure(33)
        subplot(2,2,4)
        
        hold on
        
        scatter(aLyfp*px,aIyfp/atuscal,30,green,'x','MarkerEdgeAlpha',0.15,'MarkerFaceAlpha',0.15);
        errorbar(ybin,MeanMy,StdEMy,'vertical','LineStyle','-','LineWidth',LW,'Color',green)
%         fity = fit(ybin(1:end-3)',MeanMy(1:end-3),'poly1');
%         fitplot = plot(fity,'g');

        scc = 1;
        h2=scatter(ybin,MeanMy,NumMy*scc,green,'LineWidth',LW);
        scatter(ybin,MeanMy,NumMy*scc,green,'filled','LineWidth',LW,'MarkerFaceAlpha',0.2)
        
        vlplot = vline(bins,'k:');
        set(vlplot,'LineWidth',0.5)

        hold off
        legend off
        xlabel(CL)
        ylabel('Number of proteins')
%         set(fitplot,'LineWidth',2)
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS2)
        set(gcf,'Position',pos)
        legend(h2,'Tus')
        set(gca,'XLim',[2.3,maxv])
        Ylim = get(gca,'Ylim');
        Ylim(2) = Ylim(2)/2;
        set(gca,'Ylim',Ylim)
        TitlePos = GetTitPos(Tx,Ty);t=title('D.','Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS2)
        
        % DnaN

        figure(33)
        subplot(2,2,3)
        
        hold on

        scatter(aLrfp*px,aIrfp/adnancal,30,'rx','MarkerEdgeAlpha',0.15,'MarkerFaceAlpha',0.15);
        errorbar(rbin,MeanMr,StdEMr,'vertical','LineStyle','-','LineWidth',LW,'Color','r')
%         fitr = fit(rbin(1:end-3)',MeanMr(1:end-3),'poly1');
%         fitplot = plot(fitr,'r');

        scc = 1.5;
        h3=scatter(rbin,MeanMr,NumMr*scc,'r','LineWidth',LW);
        scatter(rbin,MeanMr,NumMr*scc,'r','filled','LineWidth',LW,'MarkerFaceAlpha',0.2)
        
        vlplot = vline(bins,'k:');
        set(vlplot,'LineWidth',0.5)

        legend off
        hold off
        xlabel(CL)
        ylabel('Number of proteins')
%         set(fitplot,'LineWidth',2)
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS2)
        set(gcf,'Position',screensize3)
        legend(h3,'DnaN')
        Ylim = get(gca,'Ylim');
        Ylim(2) = Ylim(2)/3*2;
        set(gca,'Ylim',Ylim)
        set(gca,'XLim',[2.3,maxv])
        TitlePos = GetTitPos(Tx,Ty);t=title('C.','Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS2)
    end
    
    
    
%% CFP
    if c_PvsCL==1
        figure
        set(gcf, 'position', screensize2)

        % Kymo
        subplot(2,1,1)
        hold on
        vl = vline(0.9,'k--');
        set(vl,'LineWidth',LW)
        scatter(single(kLcfp),kPcfp,kIcfp/k_Intensityval(1),blue,'filled');
        % myfit=polyfit(Acfp,Bcfp,4);
        % x=0:0.01:1;
        % y=polyval(myfit,x);
        % plot(x,y,'r','LineWidth',5)
      
        hold off
        ylabel(P); xlabel(CC);
        axis([0 1 0 1])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS,'Color',[0.95 0.95 0.95])
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{1},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
        
        % Agar
        subplot(2,1,2)
        hold on
        vl = vline(4.2,'k--');
        set(vl,'LineWidth',LW)
        scatter(single(aLcfp)*px,aPcfp,aIcfp/aIntensityval(1),blue,'filled');
        myfit=polyfit(aLcfp*px,aPcfp,4);
        x=(15:0.1:45)*px;
        y=polyval(myfit,x);
        %plot(x,y,'r','LineWidth',5)
        ylabel(P); xlabel(CL);
        hold off
        axis([CLbound*px 0 1])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{2},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
    end
    
    if c_IvsP==1
        figure
        set(gcf, 'position', screensize2)
        
        % Kymo
        subplot(2,1,1)
        hold on
        scatter(kPcfp,kIcfp,'x','MarkerEdgeColor',blue);
        myfit=polyfit(kPcfp,kIcfp,4);
        x=0:00.1:1;
        y=polyval(myfit,x);
        plot(x,y,'k','LineWidth',3)
        ylabel(SI); xlabel(P);
        hold off
        axis([0 1 -0.1 35])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS,'Color',[0.95 0.95 0.95])
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{1},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
        
        % Agar
        subplot(2,1,2)
        hold on
        scatter(aPcfp,aIcfp,'x','MarkerEdgeColor',blue);
        myfit=polyfit(aPcfp,aIcfp,4);
        x=0:0.001:1;
        y=polyval(myfit,x);
        plot(x,y,'k','LineWidth',3)
        ylabel(SI); xlabel(P);
        hold off
        axis([0 1 -0.1 90])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{2},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
    end
    
    if c_NSvsP==1
        figure
        set(gcf, 'position', screensize2)
        bins = 15;
        thisedge = (0:bins)/bins;
        red = [.8,0,0];
        
        % Kymo
        subplot(2,1,1)
        [numbin,edges] = histcounts(kPcfp,20);
        X = diff(edges);
        X = cumsum(X) - X(1)/2;
        hold on
        scatter(kPcfp,kIcfp,'x','MarkerEdgeColor',blue);
        
        yyaxis right
        plot(X,numbin,'Color',red,'LineWidth',3)
        set(gca,'YColor',red)
        ylabel(NS)
        yyaxis left
        ylabel(SI); xlabel(P)
        hold off
        
        set(gca,'XLim',[0,1],'Fontname',Font,'LineWidth',AW,'FontSize',FS,'Color',[0.95 0.95 0.95])
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{1},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
        
        % Agar
        subplot(2,1,2)
        [numbin,edges] = histcounts(aPcfp,thisedge);
        X = diff(edges);
        X = cumsum(X) - X(1)/2;
        hold on
        scatter(aPcfp,aIcfp,'x','MarkerEdgeColor',blue);
        
        yyaxis right
        plot(X,numbin,'Color',red,'LineWidth',3)
        set(gca,'YColor',red)
        ylabel(NS)
        yyaxis left
        ylabel(SI); xlabel(P)
        hold off
        
        set(gca,'XLim',[0,1],'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{2},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
    end
    
    if c_Heatmap==1
        bins = 15;
        kthisedge2{1} = linspace(min(kLcfp),max(kLcfp),bins+5);
        kthisedge2{2} = (0:bins)/bins;
        athisedge2{1} = linspace(min(aLcfp*px),max(aLcfp*px),bins+5);
        athisedge2{2} = (0:bins)/bins;
        
        figure
        set(gcf, 'position', screensize2)
        colormap(hot)

        % Kymo
        kNumcfp(1,:) = kLcfp;
        kNumcfp(2,:) = kPcfp;

        subplot(2,1,1)
        Heatmap = hist3(kNumcfp','Edges',kthisedge2);
        pcolor(kthisedge2{1},(kthisedge2{2}),Heatmap');
        ylabel(P); xlabel(CC);

        grid on
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{1},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
        
        % Agar
        aNumcfp(1,:) = aLcfp*px;
        aNumcfp(2,:) = aPcfp;
        
        [bin_cfp,idx_cfp]=find(kIcfp>6000); 

        NumCFP(:,idx_cfp)=aNumcfp(:,idx_cfp);
        NumCFPnz(1,:)=nonzeros(NumCFP(1,:));
        NumCFPnz(2,:)=nonzeros(NumCFP(2,:));

        FilteredSignalCFP=size(NumCFPnz,2)/size(kIcfp,2); % Percentage of total signal

        subplot(2,1,2)
        Heatmap = hist3(NumCFPnz','Edges',athisedge2);
        pcolor(athisedge2{1},(athisedge2{2}),Heatmap');

        xlabel(CL);ylabel(P)
        
        grid on
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{2},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
    end

%% YFP
    if y_PvsCL==1
        figure
        set(gcf, 'position', screensize2)
        
        % Kymo
        subplot(2,1,1)
        hold on
        vl = vline(0.9,'k--');
        set(vl,'LineWidth',LW)
        scatter(single(kLyfp),kPyfp,kIyfp/k_Intensityval(2),green,'filled');
        % myfit=polyfit(Ayfp,Byfp,4);
        % x=0:0.01:1;
        % y=polyval(myfit,x);
        % plot(x,y,'r','LineWidth',5)
        hold off
        ylabel(P); xlabel(CC);

        axis([0 1 0 1])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS,'Color',[0.95 0.95 0.95])
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{1},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
        
        % Agar
        subplot(2,1,2)
        hold on
        vl = vline(4.2,'k--');
        set(vl,'LineWidth',LW)
        scatter(single(aLyfp*px),aPyfp,aIyfp/aIntensityval(2),green,'filled');
        myfit=polyfit(aLyfp*px,aPyfp,4);
        x=15:0.1:45;
        y=polyval(myfit,x);
        %plot(x,y,'r','LineWidth',5)
        ylabel(P); xlabel(CL);
        hold off
        axis([CLbound*px 0 1])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{2},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
    end
    
    if y_IvsP==1
        figure
        set(gcf, 'position', screensize2)
        
        % Kymo
        subplot(2,1,1)
        hold on
        scatter(kPyfp,kIyfp,'x','MarkerEdgeColor',green);
        myfit=polyfit(kPyfp,kIyfp,4);
        x=0:00.1:1;
        y=polyval(myfit,x);
        plot(x,y,'k','LineWidth',3)
        ylabel(SI); xlabel(P);
        hold off
        axis([0 1 -0.1 35])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS,'Color',[0.95 0.95 0.95])
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{1},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
        
        % Agar
        subplot(2,1,2)
        hold on
        scatter(aPyfp,aIyfp,'x','MarkerEdgeColor',green);
        myfit=polyfit(aPyfp,aIyfp,4);
        x=0:0.001:1;
        y=polyval(myfit,x);
        plot(x,y,'k','LineWidth',3)
        ylabel(SI); xlabel(P);
        hold off
        axis([0 1 -0.1 90])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{2},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
    end
    
    if y_NSvsP==1
        figure
        set(gcf, 'position', screensize2)
        bins = 15;
        thisedge = (0:bins)/bins;
        
        % Kymo
        subplot(2,1,1)
        [numbin,edges] = histcounts(kPyfp,20);
        norm = max(numbin)/35;
        X = diff(edges);
        X = cumsum(X) - X(1)/2;
        hold on
        scatter(kPyfp,kIyfp,'x','MarkerEdgeColor',green);
        yyaxis right
        plot(X,numbin,'Color',blue,'LineWidth',3)
        set(gca,'YColor',blue)
        ylabel(NS)
        yyaxis left
        ylabel(SI); xlabel(P)
        hold off
        set(gca,'XLim',[0,1],'Fontname',Font,'LineWidth',AW,'FontSize',FS,'Color',[0.95 0.95 0.95])
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{1},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
        
        % Agar
        subplot(2,1,2)
        [numbin,edges] = histcounts(aPyfp,thisedge);
        norm = max(numbin)/80;
        X = diff(edges);
        X = cumsum(X) - X(1)/2;
        hold on
        scatter(aPyfp,aIyfp,'x','MarkerEdgeColor',green);
        yyaxis right
        plot(X,numbin,'Color',blue,'LineWidth',3)
        set(gca,'YColor',blue)
        ylabel(NS)
        yyaxis left
        ylabel(SI); xlabel(P)
        hold off
        set(gca,'XLim',[0,1],'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{2},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
    end
    
    if y_Heatmap==1
        bins = 15;
        kthisedge2{1} = linspace(min(kLyfp),max(kLyfp),bins+5);
        kthisedge2{2} = (0:bins)/bins;
        athisedge2{1} = linspace(min(aLyfp*px),max(aLyfp*px),bins+5);
        athisedge2{2} = (0:bins)/bins;
        
        figure
        set(gcf, 'position', screensize2)
        colormap(hot)

        % Kymo
        kNumyfp(1,:) = kLyfp;
        kNumyfp(2,:) = kPyfp;

        subplot(2,1,1)
        Heatmap = hist3(kNumyfp','Edges',kthisedge2);
        pcolor(kthisedge2{1},(kthisedge2{2}),Heatmap');
        ylabel(P); xlabel(CC);
        grid on
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{2},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
        
        % Agar
        aNumyfp(1,:) = aLyfp*px;
        aNumyfp(2,:) = aPyfp;
        
        [bin_yfp,idx_yfp]=find(kIyfp>yfp.filterval); 

        NumYFP(:,idx_yfp)=aNumyfp(:,idx_yfp);
        NumYFPnz(1,:)=nonzeros(NumYFP(1,:));
        NumYFPnz(2,:)=nonzeros(NumYFP(2,:));

        FilteredSignalYFP=size(NumYFPnz,2)/size(kIyfp,2); % Percentage of total signal

        subplot(2,1,2)
        Heatmap = hist3(NumYFPnz','Edges',athisedge2);
        pcolor(athisedge2{1},(athisedge2{2}),Heatmap');

        ylabel(P); xlabel(CL);
        grid on
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{2},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
    end

%% RFP
    if r_PvsCL==1
        figure
        set(gcf, 'position', screensize2)
        
        % Kymo
        subplot(2,1,1)
        hold on
        vl = vline(0.9,'k--');
        set(vl,'LineWidth',LW)
        scatter(single(kLrfp),kPrfp,kIrfp/k_Intensityval(3),'r','filled');
        % myfit=polyfit(Arfp,Brfp,4);
        % x=0:0.01:1;
        % y=polyval(myfit,x);
        % plot(x,y,'r','LineWidth',5)
        hold off
        ylabel(P); xlabel(CC);
        axis([0 1 0 1])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS,'Color',[0.95 0.95 0.95])
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{1},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
        
        % Agar
        subplot(2,1,2)
        hold on
        vl = vline(4.2,'k--');
        set(vl,'LineWidth',LW)
        scatter(single(aLrfp*px),aPrfp,aIrfp/aIntensityval(3),'r','filled');
        myfit=polyfit(aLrfp*px,aPrfp,4);
        x=15:0.1:45;
        y=polyval(myfit,x);
        %plot(x,y,'r','LineWidth',5)
         ylabel(P); xlabel(CL);
        hold off
        axis([CLbound*px 0 1])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{2},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
    end
    
    if r_IvsP==1
        figure
        set(gcf, 'position', screensize2)
        
        % Kymo
        subplot(2,1,1)
        hold on
        scatter(kPrfp,kIrfp,'x','MarkerEdgeColor','r');
        myfit=polyfit(kPrfp,kIrfp,4);
        x=0:00.1:1;
        y=polyval(myfit,x);
        plot(x,y,'k','LineWidth',3)
        ylabel(SI); xlabel(P);
        hold off
        axis([0 1 -0.1 35])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS,'Color',[0.95 0.95 0.95])
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{1},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
        
        % Agar
        subplot(2,1,2)
        hold on
        scatter(aPrfp,aIrfp,'x','MarkerEdgeColor','r');
        myfit=polyfit(aPrfp,aIrfp,4);
        x=0:0.001:1;
        y=polyval(myfit,x);
        plot(x,y,'k','LineWidth',3)
        ylabel(SI); xlabel(P);
        hold off
        axis([0 1 -0.1 90])
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{2},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
    end
    
    if r_NSvsP==1
        figure
        set(gcf, 'position', screensize2)
        bins = 15;
        thisedge = (0:bins)/bins;
        
        % Kymo
        subplot(2,1,1)
        [numbin,edges] = histcounts(kPrfp,20);
        norm = max(numbin)/35;
        X = diff(edges);
        X = cumsum(X) - X(1)/2;
        hold on
        scatter(kPrfp,kIrfp,'x','MarkerEdgeColor','r');
        
        yyaxis right
        plot(X,numbin,'Color',blue,'LineWidth',3)
        set(gca,'YColor',blue)
        ylabel(NS)
        yyaxis left
        ylabel(SI); xlabel(P)
        hold off
        
        set(gca,'XLim',[0,1],'Fontname',Font,'LineWidth',AW,'FontSize',FS,'Color',[0.95 0.95 0.95])
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{1},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
        
        % Agar
        subplot(2,1,2)
        [numbin,edges] = histcounts(aPrfp,thisedge);
        norm = max(numbin)/80;
        X = diff(edges);
        X = cumsum(X) - X(1)/2;
        hold on
        scatter(aPrfp,aIrfp,'x','MarkerEdgeColor','r');
        
        yyaxis right
        plot(X,numbin,'Color',blue,'LineWidth',3)
        set(gca,'YColor',blue)
        ylabel(NS)
        yyaxis left
        ylabel(SI); xlabel(P)
        hold off
        
        set(gca,'XLim',[0,1],'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{2},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
    end
    
    if r_Heatmap==1
        bins = 15;
        kthisedge2{1} = linspace(min(kLrfp),max(kLrfp),bins+5);
        kthisedge2{2} = (0:bins)/bins;
        athisedge2{1} = linspace(min(aLrfp*px),max(aLrfp*px),bins+5);
        athisedge2{2} = (0:bins)/bins;
        
        figure
        set(gcf, 'position', screensize2)
        colormap(hot)

        % Kymo
        kNumrfp(1,:) = kLrfp;
        kNumrfp(2,:) = kPrfp;

        subplot(2,1,1)
        Heatmap = hist3(kNumrfp','Edges',kthisedge2);
        pcolor(kthisedge2{1},(kthisedge2{2}),Heatmap');
        ylabel(P); xlabel(CC);
        grid on
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{1},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
        
        % Agar
        aNumrfp(1,:) = aLrfp*px;
        aNumrfp(2,:) = aPrfp;

        [bin_rfp,idx_rfp]=find(kIrfp>10000); 

        NumRFP(:,idx_rfp)=aNumrfp(:,idx_rfp);
        NumRFPnz(1,:)=nonzeros(NumRFP(1,:));
        NumRFPnz(2,:)=nonzeros(NumRFP(2,:));

        FilteredSignalRFP=size(NumRFPnz,2)/size(kIrfp,2);

        subplot(2,1,2)
        Heatmap = hist3(NumRFPnz','Edges',athisedge2);
        h = pcolor(athisedge2{1},(athisedge2{2}),Heatmap');

        ylabel(P); xlabel(CL);
        grid on
        set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
        TitlePos = GetTitPos(Tx,Ty);t=title(Title{2},'Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
    end


    
%%

if all_Heatmap==1
    bins = 20;
    
    figure
    set(gcf, 'position', screensize2)

    %% CFP
    kthisedge2{1} = linspace(min(kLcfp),max(kLcfp),bins+5);
    kthisedge2{2} = (0:bins)/bins;
    athisedge2{1} = linspace(min(aLcfp*px),max(aLcfp*px),bins+5);
    athisedge2{2} = (0:bins)/bins;

    % Kymo
    kNumcfp(1,:) = kLcfp;
    kNumcfp(2,:) = kPcfp;

    ax1=subplot(3,2,1);
    Heatmap = hist3(kNumcfp','Edges',kthisedge2);
    pcolor(kthisedge2{1},(kthisedge2{2}),Heatmap');

    grid on
    set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
    colormap(ax1,winter)
    colorbar('FontName',Font,'LineWidth',AW,'FontSize',FS)
    TitlePos = GetTitPos(Tx,Ty);t=title('A.','Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)

    % Agar
    aNumcfp(1,:) = aLcfp*px;
    aNumcfp(2,:) = aPcfp;

    [bin_cfp,idx_cfp]=find(aIcfp>6000); 

    NumCFP(:,idx_cfp)=aNumcfp(:,idx_cfp);
    NumCFPnz(1,:)=nonzeros(NumCFP(1,:));
    NumCFPnz(2,:)=nonzeros(NumCFP(2,:));

    FilteredSignalCFP=size(NumCFPnz,2)/size(aIcfp,2); % Percentage of total signal

    ax2=subplot(3,2,2);
    Heatmap = hist3(NumCFPnz','Edges',athisedge2);
    pcolor(athisedge2{1},(athisedge2{2}),Heatmap');

    grid on
    set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
    colormap(ax2,winter)
    TitlePos = GetTitPos(Tx,Ty);t=title('B.','Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
    
    
    %% YFP
    kthisedge2{1} = linspace(min(kLyfp),max(kLyfp),bins+5);
    kthisedge2{2} = (0:bins)/bins;
    athisedge2{1} = linspace(min(aLyfp*px),max(aLyfp*px),bins+5);
    athisedge2{2} = (0:bins)/bins;

    % Kymo
    kNumyfp(1,:) = kLyfp;
    kNumyfp(2,:) = kPyfp;

    ax3=subplot(3,2,3);
    Heatmap = hist3(kNumyfp','Edges',kthisedge2);
    pcolor(kthisedge2{1},(kthisedge2{2}),Heatmap');
    ylabel(P);
    grid on
    set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
    colormap(ax3,summer)
    colorbar('FontName',Font,'LineWidth',AW,'FontSize',FS)
    TitlePos = GetTitPos(Tx,Ty);t=title('C.','Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)

    % Agar
    aNumyfp(1,:) = aLyfp*px;
    aNumyfp(2,:) = aPyfp;

    [bin_yfp,idx_yfp]=find(aIyfp>0.5*yfp.filterval); 

    NumYFP(:,idx_yfp)=aNumyfp(:,idx_yfp);
    NumYFPnz(1,:)=nonzeros(NumYFP(1,:));
    NumYFPnz(2,:)=nonzeros(NumYFP(2,:));

    FilteredSignalYFP=size(NumYFPnz,2)/size(aIyfp,2); % Percentage of total signal

    ax4=subplot(3,2,4);
    Heatmap = hist3(NumYFPnz','Edges',athisedge2);
    pcolor(athisedge2{1},(athisedge2{2}),Heatmap');

    grid on
    set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
    colormap(ax4,summer)
    TitlePos = GetTitPos(Tx,Ty);t=title('D.','Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
    
    %% RFP
    
    kthisedge2{1} = linspace(min(kLrfp),max(kLrfp),bins+5);
    kthisedge2{2} = (0:bins)/bins;
    athisedge2{1} = linspace(min(aLrfp*px),max(aLrfp*px),bins+5);
    athisedge2{2} = (0:bins)/bins;

    % Kymo
    kNumrfp(1,:) = kLrfp;
    kNumrfp(2,:) = kPrfp;

    ax5=subplot(3,2,5);
    Heatmap = hist3(kNumrfp','Edges',kthisedge2);
    pcolor(kthisedge2{1},(kthisedge2{2}),Heatmap');
    xlabel(CC);
    grid on
    set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
    colormap(ax5,hot)
    colorbar('FontName',Font,'LineWidth',AW,'FontSize',FS)
    TitlePos = GetTitPos(Tx,Ty);t=title('E.','Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)

    % Agar
    aNumrfp(1,:) = aLrfp*px;
    aNumrfp(2,:) = aPrfp;

    [bin_rfp,idx_rfp]=find(aIrfp>10000); 

    NumRFP(:,idx_rfp)=aNumrfp(:,idx_rfp);
    NumRFPnz(1,:)=nonzeros(NumRFP(1,:));
    NumRFPnz(2,:)=nonzeros(NumRFP(2,:));

    FilteredSignalRFP=size(NumRFPnz,2)/size(aIrfp,2);

    ax6=subplot(3,2,6);
    Heatmap = hist3(NumRFPnz','Edges',athisedge2);
    h = pcolor(athisedge2{1},(athisedge2{2}),Heatmap');

    xlabel(CL);
    grid on
    set(gca,'Fontname',Font,'LineWidth',AW,'FontSize',FS)
    colormap(ax6,hot)
    TitlePos = GetTitPos(Tx,Ty);t=title('F.','Position',TitlePos);set(t,'FontName',TitleFont,'FontSize',TitleFS)
end

