%FRAP processing
%Tristan Ursell
%August 2011

function frap_processing

[fname,pname,~]=uigetfile({'*.tif','TIF Image File'});

if fname==0
    return
end

N=length(imfinfo([pname fname]));

h1=figure;
imagesc(imread([pname fname],1))
title('Choose FRAPed region.')
axis equal
axis tight
colormap(hot)
rect1=round(getrect);
close(h1)

h2=figure;
imagesc(imread([pname fname],1))
title('Choose non-FRAPed region.')
axis equal
axis tight
colormap(hot)
rect2=round(getrect);
close(h2)

prompt={'Enter time step between frames(s): '};
name='Temporal Calibration';
numlines=1;
defaultanswer={''};
answer=inputdlg(prompt,name,numlines,defaultanswer);

if and(~isempty(answer),~isempty(str2num(answer{1})))
    dt=str2num(answer{1});
else
    dt=1;
end
timevec=0:dt:(N-1)*dt;
    
frp_mean=zeros(1,N);
frp_std=zeros(1,N);
nfrp_mean=zeros(1,N);
nfrp_std=zeros(1,N);
for i=1:N
    I0=double(imread([pname fname],i));
    frp_temp=I0(rect1(2):rect1(2)+rect1(4),rect1(1):rect1(1)+rect1(3));
    nfrp_temp=I0(rect2(2):rect2(2)+rect2(4),rect2(1):rect2(1)+rect2(3));
    
    frp_mean(i)=mean(frp_temp(:));
    frp_std(i)=std(frp_temp(:));
    nfrp_mean(i)=mean(nfrp_temp(:));
    nfrp_std(i)=std(nfrp_temp(:));
end

figure;
hold on
plot(timevec,frp_mean,'Linewidth',2,'color',[0.9,0,0])
plot(timevec,frp_mean-frp_std,'--','Linewidth',2,'color',[0.9,0,0])
plot(timevec,frp_mean+frp_std,'--','Linewidth',2,'color',[0.9,0,0])

plot(timevec,nfrp_mean,'Linewidth',2,'color',[0,0.9,0])
plot(timevec,nfrp_mean-nfrp_std,'--','Linewidth',2,'color',[0,0.9,0])
plot(timevec,nfrp_mean+nfrp_std,'--','Linewidth',2,'color',[0,0.9,0])

xlabel('Time(s)')
ylabel('Intensity')
box on

figure;
hold on
plot([0,N*dt],[1,1],'k--')
plot(timevec,frp_mean./nfrp_mean,'Linewidth',2,'color',[0.9,0,0])
xlabel('Time(s)')
ylabel('Relative Recovery')
xlim([0,N*dt])
ymax=(ceil(10*max(frp_mean./nfrp_mean))+1)/10;
ylim([0,ymax+0.001])
box on

%Fitting relative recovoery
dataq = questdlg('Would you like to fit the relative recovery data?', ...
    'Recovery Fitting', ...
    'Yes', 'No', 'Yes');

if strcmp('Yes',dataq)
    relrecov=frp_mean./nfrp_mean;
    frapstart=min(find(diff(relrecov)==min(diff(relrecov))))+1;
    Xdata=timevec(frapstart:end);
    Ydata=relrecov(frapstart:end);
    
    g=@(x) sum((Ydata-(1-x(1)*exp(-Xdata/x(2)))).^2);
    x0=[0.5,2*(max(Ydata)-min(Ydata))/(max(Xdata)-min(Xdata))];
    outp=fminsearch(g,x0);
    
    plot(Xdata,1-outp(1)*exp(-Xdata/outp(2)),'k','Linewidth',2)
    title(['Time constant(s): ' num2str(outp(2))])
end




