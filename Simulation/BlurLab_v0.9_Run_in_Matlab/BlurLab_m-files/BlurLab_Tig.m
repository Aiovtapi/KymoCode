function varargout = BlurLab_Tig(varargin)
%
% For instructions, launch BlurLab and select the 'Help' button.
%
% The following m-files should be included with 'BlurLab_Tig.m' and
% 'BlurLab_Tig.fig':
%
% custom_noise.m
% helppopup.m
% histupdate_Callback.m
% BlurLab_text.m
% loaddefaults_Callback.m
% loadpars_Callback.m
% matoverlay.m
% matoverlay_sub.m
% memupdate_Callback.m
% modone.m
% parametersave_Callback.m
% poissimnoise.m
% pushbutton_state_Callback.m
% random_box.m
%
% There should also be:
%
% BlurLab_Help.txt
% BlurLab_License.txt
% PSF_100X_na1.4_n1.515_lam532_dx80_dz10.tif
% MCRInstaller.exe
% BlurLab_icon.jpg
%
% Optional files for creating input text files with Matlab are included in
% folder 'M-File_Models'.
%

% Last Modified by GUIDE v2.5 01-Jul-2016 13:36:21

% Skip frame option?
% only frap certain Z thickness

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @BlurLab_Tig_OpeningFcn, ...
    'gui_OutputFcn',  @BlurLab_Tig_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before BlurLab_Tig is made visible.
function BlurLab_Tig_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;

%set checks
global loadch
global psfch
global simch
global frapch
global ncols
global noisech
global cd_write
loadch=0;
psfch=0;
simch=0;
noisech=0;
frapch=0;
ncols=5;
cd_write=0;

set(handles.loaddefaults,'Visible','off')
set(handles.memupdate,'Visible','off')
set(handles.histupdate,'Visible','off')
set(handles.text_boxcar,'Visible','off')

set(handles.text49,'Enable','off')
set(handles.slider1,'Enable','off')
set(handles.edit_znum,'Enable','off')
set(handles.text_slidemax,'Enable','off')
set(handles.text_slidemin,'Enable','off')

set(handles.uipanel_output,'Visible','off')
set(handles.uipanel_size,'Visible','off')
set(handles.uipanel_frames,'Visible','off')
set(handles.uipanel_zstacking,'Visible','off')
set(handles.uipanel_randinput,'Visible','off')
set(handles.uipanel_frap,'Visible','off')
set(handles.uipanel_noise,'Visible','off')
set(handles.uipanel_subnoise,'Visible','off')

%clear temp image file
if ~isempty(dir('BlurLab_temp.tif'))
    delete BlurLab_temp.tif
end

%format axes
axes(handles.axes1);
axis equal
xlabel('x (um)')
ylabel('y (um)')
box on

memupdate_Callback(hObject, eventdata, handles)
loaddefaults_Callback(hObject, eventdata, handles)
checkbox_histogram_Callback(hObject, eventdata, handles)
popupmenu1_Callback(hObject, eventdata, handles)
guidata(hObject,handles)

% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)

% --- Outputs from this function are returned to the command line.
function varargout = BlurLab_Tig_OutputFcn(hObject, eventdata, handles)
varargout{1} = handles.output;

% --- Executes on button press in pushbutton_sim.
function pushbutton_sim_Callback(hObject, eventdata, handles)
global basename
global Npsf
global file2
global loadch
global psfch
global frapch
global dz
global simch
global xpos
global ypos
global Ikeep
global textput
global sx2
global sy2
global cd1
global pvals
global intvals
global ncols
global cd_write

memupdate_Callback(hObject, eventdata, handles)

writeq=get(handles.checkbox_write,'Value');
stochq=get(handles.checkbox_stochastic,'Value');
flipq=get(handles.checkbox_flip,'Value');
frapq=get(handles.checkbox_frap,'Value');
zstackq=get(handles.checkbox_simZ,'Value');
startim=str2num(get(handles.edit_startim,'String'));
endim=str2num(get(handles.edit_endim,'String'));
basename=get(handles.edit_fname,'String');
dz=str2num(get(handles.edit_psfdz,'String'));
outX=str2num(get(handles.edit_outX,'String'));
outY=str2num(get(handles.edit_outY,'String'));
npoints=str2num(get(handles.edit_points,'String'));
imcalb=str2num(get(handles.edit_impxlsize,'String'));
sfc=str2num(get(handles.edit_sfc,'String'));
frapframe=str2num(get(handles.edit_frapframe,'String'));
boxcarq=get(handles.checkbox_boxcar,'Value');
box_frames=str2num(get(handles.edit_boxcar,'String'));

%make sure files are loaded
if psfch~=1
    warndlg('Load PSF file to begin simulating.','Cannot simulate.')
    return
elseif loadch~=1
    warndlg('Load input file or select Random Input to begin simulating.','Cannot simulate.')
    return
end

%set image minima and maxima
Imins=2^16-1;
Imaxs=0;

if and(frapch==0,frapq==1)
    set(handles.checkbox_frap,'Value',0)
    set(handles.popupmenu1,'Value',6)
    guidata(hObject,handles)
    warndlg('To use FRAP, first select a FRAP region.','FRAP Invalid')
    return
elseif frapq==1
    global xrectmin
    global xrectmax
    global yrectmin
    global yrectmax
end

%make sure if stochastic bleach, 6 columns
if and(stochq==1,ncols==5)
    set(handles.checkbox_stochastic,'Value',0)
    stochq=0;
    warndlg('To use stochastic bleaching, labels must be included with the input file.  The simulation will run without stochastic bleaching.','Stochastic Bleaching Error')
end

%get Z-stack parameters
minz=str2num(get(handles.edit_minz,'String'));
maxz=str2num(get(handles.edit_maxz,'String'));
dzstack=str2num(get(handles.edit_dzstack,'String'));
if and(zstackq==1,and(~isempty(dzstack),and(~isempty(minz),~isempty(minz))))
    zslice=linspace(minz,maxz,dzstack);
else
    zslice=str2num(get(handles.edit_zslice,'String'));
end

%choose step size for simulation loop w/ Boxcar averaging
if boxcarq==1
    loopstep=box_frames;
    set(handles.text_boxcar,'Visible','on')
    set(handles.text_boxcar,'String','Boxcar averaging: ')
    set(handles.slider1,'Min',ceil(startim/loopstep))
    set(handles.slider1,'Max',ceil(endim/loopstep))
    set(handles.text_slidemin,'String',num2str(ceil(startim/loopstep)))
    set(handles.text_slidemax,'String',num2str(ceil(endim/loopstep)))
    set(handles.slider1,'Value',ceil(startim/loopstep))
    set(handles.edit_slidenum,'String',num2str(ceil(startim/loopstep)))
else
    loopstep=1;
end

%total number of simulated images
numims=(endim-startim+1)*length(zslice)/loopstep;

%construct list matching psf Z values
Zvals=dz*(-floor(Npsf/2):1:floor(Npsf/2));

%setup holding matrix
if get(handles.checkbox_ram,'Value')
    try
        if flipq==0
            Ikeep=zeros(outY,outX,numims);
        else
            Ikeep=zeros(outX,outY,numims);
        end
    catch
        if flipq==0
            Ikeep=zeros(outY,outX);
        else
            Ikeep=zeros(outX,outY);
        end
        set(handles.checkbox_disk,'Value',1)
        set(handles.checkbox_ram,'Value',0)
        warndlg('Your computer does not support matrices this large in RAM.  BlurLab will use disk space.','Not Enough Memory');
        uiwait
    end
end

memq=get(handles.checkbox_ram,'Value');

%clear temp image file
if ~isempty(dir('BlurLab_temp.tif'))
    delete BlurLab_temp.tif
end

%adjust textput with offsets
textput_off=textput;
textput_off(:,1)=textput(:,1)+str2num(get(handles.edit_offsetx,'String'));
textput_off(:,2)=textput(:,2)+str2num(get(handles.edit_offsety,'String'));
nlines=length(textput(:,1));

%assign cdout directory for image writing
if cd_write==0
    cdout=cd1;
else
    cdout=cd_write;
end

%get output bitdepth
if get(handles.checkbox_8bit,'Value')==1
    bitdepth=8;
else
    bitdepth=16;
end

%check for redundant file output names and create unique, sequential names
if writeq==1
    l1=dir([cdout basename]);
    if ~isempty(l1)
        dotp=max(find(basename=='.'));
        numpos=max(findstr(basename,'_n'));
        if ~isempty(numpos)
            fnum=1+str2num(basename(numpos+2:dotp-1));
            basename=[basename(1:numpos-1) '_n' num2str(fnum) basename(dotp:end)];
        else
            basename=[basename(1:dotp-1) '_n1' basename(dotp:end)];
        end
        set(handles.edit_fname,'String',basename)
    end
end

%write output parameters file
if get(handles.checkbox_savepars,'Value')
    dotp=max(find(basename=='.'));
    parname=[basename(1:dotp) 'txt'];
    parametersave_Callback(hObject, eventdata, handles, cdout, parname)
end

%load psf file
Ipsf=zeros(sy2,sx2,Npsf,'single');
h3=waitbar(0,'Loading PSF ...');
for i=1:Npsf
    Ipsf(:,:,i)=imread(file2,i);
    waitbar(i/Npsf,h3)
end
close(h3);
pause(0.01)

%normalize the PSF
if get(handles.checkbox_psfnorm,'Value')==1
    minpsf=min(min(min(Ipsf)));
    maxpsf=max(max(max(Ipsf)));
    Ipsf=single(Ipsf/maxpsf);
end

%get pixel calibration
calb=str2num(get(handles.edit_impxlsize,'String'));

%initiate vectors for tracking who is bleached/frapped
if ncols==6
    bleacherslabel=zeros(1,max(textput_off(:,6)));
    fraplistlabel=zeros(1,max(textput_off(:,6)));
end

%memory update
memupdate_Callback(hObject, eventdata, handles)

%mark simulation
simch=1;

%initiate frame dlist tracking 
herelength=zeros(1,max(textput(:,5)));

%start simulation
q=0;
for i=startim:loopstep:endim
    clear dlist Xlist Ylist Zlist Ilist Llist framelist
    %list of frames for current image
    framelist=i:(i+loopstep-1);
        
    %find data for current frame with intensity greater than zero
    dtemp=zeros(1,nlines);
    nl=0;
    for k=framelist
        dtemp0=find(and(textput_off(:,5)==k,textput_off(:,4)>0));
        dtemp(1+nl:nl+length(dtemp0))=dtemp0;
        nl=nl+length(dtemp0);
        herelength(k)=nl;
    end
    dlist=dtemp(1:nl);
    
    %list of matrix positions
    Xlist=textput_off(dlist,1);
    Ylist=textput_off(dlist,2);
    Zlist=textput_off(dlist,3);
    
    %get label list
    if ncols==6
        Llist=textput_off(dlist,6);
    end
    
    %handle TIRF
    if get(handles.checkbox_tirf,'Value')
        lam_tirf=str2num(get(handles.edit_tirfdecay,'String'));
        zoff_tirf=str2num(get(handles.edit_tirfZ,'String'));
        Ilist=exp(-(Zlist-zoff_tirf)/lam_tirf).*textput_off(dlist,4);
    else
        Ilist=textput_off(dlist,4);
    end
    
    %handle FRAPing
    if and(sum(framelist==frapframe)>0,and(ncols==6,get(handles.checkbox_frap,'Value')))
        %get percentage
        per=str2num(get(handles.edit_frap,'String'))/100;
        
        %find the labels of the points frapframe
        dfrap=find(textput_off(:,5)==frapframe);
        Lfrap=textput_off(dfrap,6);
        Xfrap=textput_off(dfrap,1);
        Yfrap=textput_off(dfrap,2);
        
        if get(handles.checkbox_flip,'Value')==0
            rectlist=(Xfrap>=xrectmin).*(Xfrap<=xrectmax).*...
                (Yfrap>=yrectmin).*(Yfrap<=yrectmax);
        else
            rectlist=(Yfrap>=xrectmin).*(Yfrap<=xrectmax).*...
                (Xfrap>=yrectmin).*(Xfrap<=yrectmax);
        end
        fraplist=(rand(length(Xfrap),1)<=per).*rectlist;
        fraplistlabel(Lfrap(fraplist==1))=1;
    end
    
    %go through Z-slices (even if there's just one)
    for p=1:length(zslice)
        q=q+1;
        %check for stop
        global stopch
        if stopch==1
            stopch=0;
            return
        end
        
        %initialize image
        Iout=zeros(outY,outX);
        
        %stochastic photobleaching
        if and(get(handles.checkbox_bleach,'Value'),and(ncols==6,get(handles.checkbox_stochastic,'Value')))
            tau=str2num(get(handles.edit_tau,'String'));
            pout=1-exp(-1/tau); %since the frames are combined, each
            %emitter has 'loopstep' chances to bleach, which is the
            %same as setting the exponent to '-loopstep/tau' for a
            %single averaged frame
            bleachers=pout>rand(1,length(Llist));
            bleacherslabel(Llist(bleachers))=1;
        end
        
        %interpolate Z-PSF values
        %Zassign(1,:)=lower plane number
        %Zassign(2,:)=lower plane weighting
        %Zassign(3,:)=upper plane number
        %Zassign(4,:)=upper plane weighting
        Zassign=zeros(4,length(dlist));
        
        inbnds=and((Zlist-zslice(p))>=min(Zvals),(Zlist-zslice(p))<=max(Zvals));
        for k=1:length(dlist)
            if inbnds(k)
                temp1=abs(Zvals-(Zlist(k)-zslice(p)))<dz;
                fnd=find(temp1);
                if sum(temp1)==1
                    Zassign(1,k)=fnd;
                    Zassign(2,k)=1/2;
                    Zassign(3,k)=fnd;
                    Zassign(4,k)=1/2;
                elseif sum(temp1)==2
                    Zassign(1,k)=min(fnd);
                    Zassign(2,k)=1/dz*(Zvals(max(fnd))-(Zlist(k)-zslice(p)));
                    Zassign(3,k)=max(fnd);
                    Zassign(4,k)=1/dz*((Zlist(k)-zslice(p))-Zvals(min(fnd)));
                else
                    disp(['Z-slice alignment error at row ' num2str(i) '.'])
                end
            end
        end
        
        %perform the convolution
        frct=0;
        for j=1:length(dlist)
            if Zassign(1,j)~=0
                psf_temp=Ipsf(:,:,Zassign(1,j))*Zassign(2,j)+Ipsf(:,:,Zassign(3,j))*Zassign(4,j);
                if ncols==6
                    if and(bleacherslabel(Llist(j))==0,fraplistlabel(Llist(j))==0)
                        Iout=Iout+Ilist(j)*matoverlay_sub(Iout,psf_temp,Xlist(j)/calb,Ylist(j)/calb);
                    end
                else
                    Iout=Iout+Ilist(j)*matoverlay_sub(Iout,psf_temp,Xlist(j)/calb,Ylist(j)/calb);
                end
            end
            if sum(j==herelength)
                frct=frct+1;
                set(handles.text_boxcar,'String',['Boxcar averaging:  ' num2str(frct) ' / ' num2str(loopstep)])
                pause(0.001)
            end
        end
        
        %perform average across loopsteps
        Iout=Iout/loopstep;
        
        %mean field photo-bleaching
        if and(~get(handles.checkbox_stochastic,'Value'),get(handles.checkbox_bleach,'Value'))
            tau=str2num(get(handles.edit_tau,'String'));
            Iout=Iout*exp(-(i-startim)/tau);
        end
        
        %scale intensity
        Iout=Iout*str2num(get(handles.edit_sfc,'String'));
        
        %add Poisson/Shot noise
        if get(handles.checkbox_poiss,'Value')
            lam=str2num(get(handles.edit_poiss,'String'));
            basal=str2num(get(handles.edit_basal,'String'));
            Iout=poissimnoise(uint16(Iout),lam,basal);
        end
        
        %add custom noise spectrum
        if get(handles.checkbox_custom,'Value')
            T1=custom_noise(pvals,1,outX*outY);
            Iout=Iout+reshape(intvals(T1),outY,outX);;
        end
        
        %add Gaussian noise
        if get(handles.checkbox_noise,'Value')
            mu_gauss=str2num(get(handles.edit_mean,'String'));
            std_gauss=str2num(get(handles.edit_std,'String'));
            Iout=Iout+abs(normrnd(mu_gauss,std_gauss,outY,outX));
        end
        
        %redefine image for flipping
        if flipq==1
            Itemp=Iout';
            clear Iout
            Iout=Itemp;
        end
        
        %keep images
        Imins=min([Imins,min(min(Iout))]);
        Imaxs=max([Imaxs,max(max(Iout))]);
        if memq==1
            Ikeep(:,:,q)=Iout;
        else
            er1=0;
            while er1==0
                try
                    imwrite(uint16(Iout),[cd1 'BlurLab_temp.tif'],'Compression','none','Writemode','append');
                    er1=1;
                catch
                    pause(rand)
                end
            end
        end
        
        %write output images
        if writeq==1
            er1=0;
            while er1==0
                try
                    if bitdepth==8
                        imwrite(uint8(Iout),[cdout basename],'Compression','none','Writemode','append');
                    else
                        imwrite(uint16(Iout),[cdout basename],'Compression','none','Writemode','append');
                    end
                    er1=1;
                catch
                    pause(rand)
                end
            end
        end
        
        %make output plot
        if flipq==0
            xpos=0:imcalb:(outX-1)*imcalb;
            ypos=0:imcalb:(outY-1)*imcalb;
        else
            ypos=0:imcalb:(outX-1)*imcalb;
            xpos=0:imcalb:(outY-1)*imcalb;
        end
        
        axes(handles.axes1);
        if ~get(handles.checkbox_fix,'Value')
            imagesc(xpos/1000,ypos/1000,Iout)
        else
            intlow=str2num(get(handles.edit_lowfix,'String'));
            inthigh=str2num(get(handles.edit_highfix,'String'));
            imagesc(xpos/1000,ypos/1000,Iout,[intlow,inthigh])
        end
        axis equal
        set(gca,'color', [0.1 0.1 0.1])
        xlabel('x (um)')
        ylabel('y (um)')
        box on
        
        %update histogram
        if get(handles.checkbox_histogram,'Value')
            histupdate_Callback(hObject, eventdata, handles, Iout)
        end
        
        %update slider position
        set(handles.edit_slidenum,'String',num2str(ceil(i/loopstep)))
        set(handles.slider1,'Value',ceil(i/loopstep))
        
        set(handles.text_intmin,'String',num2str(1/100*round(100*min(min(Iout)))));
        set(handles.text_intmax,'String',num2str(1/100*round(100*max(max(Iout)))));
        set(handles.edit_znum,'String',num2str(p));
        
        memupdate_Callback(hObject, eventdata, handles)
        guidata(hObject,handles)
        pause(0.01)
    end
end

set(handles.text_boxcar,'Visible','off')

if str2num(get(handles.edit_sfc,'String'))~=1
    menutemp=get(handles.popupmenu1,'Value');
    set(handles.uipanel_output,'Visible','off')
    set(handles.uipanel_size,'Visible','off')
    set(handles.uipanel_frames,'Visible','off')
    set(handles.uipanel_zstacking,'Visible','off')
    set(handles.uipanel_randinput,'Visible','off')
    set(handles.uipanel_frap,'Visible','off')
    set(handles.uipanel_noise,'Visible','on')
    set(handles.uipanel_subnoise,'Visible','on')
    for qq=1:7
        set(handles.edit_sfc,'Visible','off')
        pause(0.25)
        set(handles.edit_sfc,'Visible','on')
        pause(0.125)
    end
    set(handles.uipanel_noise,'Visible','off')
    set(handles.uipanel_subnoise,'Visible','off')
    switch menutemp
        case 1
            set(handles.uipanel_frames,'Visible','on')
        case 2
            set(handles.uipanel_size,'Visible','on')
        case 3
            set(handles.uipanel_noise,'Visible','on')
            set(handles.uipanel_subnoise,'Visible','on')
        case 4
            set(handles.uipanel_zstacking,'Visible','on')
        case 5
            set(handles.uipanel_output,'Visible','on')
        case 6
            set(handles.uipanel_frap,'Visible','on')
        case 7
            set(handles.uipanel_randinput,'Visible','on')
    end
end

% --- Executes on button press in pushbutton_load.
function pushbutton_load_Callback(hObject, eventdata, handles)
global basename
global Nfile
global sx
global sy
global loadch
global textput
global cd1
global ncols

memupdate_Callback(hObject, eventdata, handles)

[fullname,cd1,filter1]=uigetfile({'*.txt','Text Input (*.txt)';'*.tif','Image Stack (*.tif)'},'Pick a coordinate input file.');

if fullname==0
    return
end

if cd1==0
    cd1=char([]);
end

loadch=1;
[fn]=max(find(fullname=='.'));
fname=fullname(1:fn-1);
file1=[cd1 fullname];

%get offset
offx=str2num(get(handles.edit_offsetx,'String'));
offy=str2num(get(handles.edit_offsety,'String'));

%get pixel calibration
pxlconv=str2num(get(handles.edit_impxlsize,'String'));

%If it's a TIF file, convert to a text file
if filter1==2
    prompt={'Z slices per frame:','Z spacing between frames (nm):'};
    name='Frame Info';
    numlines=1;
    defaultanswer={'1','100'};
    
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    
    if isempty(answer)
        return
    end
    
    nframes=str2num(char(answer(1)));
    dzframe=str2num(char(answer(2)));
    
    if or(nframes<1,dzframe<1)
        warndlg('Z slices per frame and frame spacing must be positive values.','Bad Parameters')
        return
    end
    
    fileinfo=imfinfo(file1);
    Nfile=floor(length(fileinfo)/nframes);
    
    [sy, sx]=size(imread(file1,1));
    sizestring=[num2str(sy) ' x ' num2str(sx)];
    
    %create output text file for image spots
    fclose('all');
    txtname=['text_convert_' fname '_' date '.txt'];
    fid = fopen([cd1 txtname],'w+');
    
    %*****TEXT OUT******
    wbar= waitbar(0,'Converting images to text file ...');
    ff=0;
    for i=1:length(fileinfo)
        if modone(i,nframes)==1
            ff=ff+1;
            zcount=0;
        end
        
        %write output info to text file
        Ifind=imread(file1,i);
        
        %get xy points and intensity
        [Ys,Xs]=find(Ifind>0);
        [list]=find(Ifind>0);
        Ints=double(Ifind(list));
        
        %compute z position
        Zs=ones(length(Ints),1)*dzframe*zcount;
        
        %get frame number
        ffs=ones(length(Ints),1)*ff;
        
        outdata=[Xs*pxlconv,Ys*pxlconv,Zs,Ints,ffs];
        
        %fprintf(fid,['X(nm) Y(nm) Z(nm) Int Frame\n']);
        fprintf(fid,'%6.1f %6.1f %6.1f %6.1f %6f\n',outdata');
        
        zcount=zcount+1;
        waitbar(i/Nfile,wbar)
    end
    fclose(fid);
    close(wbar);
    
    %change file1 to the text file
    file1=[cd1 txtname];
    %load text file
    h1=msgbox('Loading text file ...','Text Load','warn');
    
    %textput=load(file1);
    textput=importdata(file1);
    %fid1=fopen(file1);
    %textput=cell2mat(textscan(fid1,'%f32'));
    %fclose(fid1);
    
    close(h1);
else
    %figure out file calibration
    button1 = questdlg('In which units are the source positions?', ...
        'Text Input Units', ...
        'nanometers', 'microns', 'nanometers');
end

%load text file
h1=warndlg('Loading text file ...','Text Load');
%texttemp=load(file1);
texttemp=importdata(file1);

%fid1=fopen(file1);
%texttemp=cell2mat(textscan(fid1,'%f32'));
%fclose(fid1);

%make sure file has correct number of columsn
if or(size(texttemp,2)==5,size(texttemp,2)==6)
    if size(texttemp,2)==5
        ncols=5;
    else
        ncols=6;
    end
else
    warndlg('Text input must be formatted with columns: X,Y, Z, Intensity, Frame, Label.','Incorrect Formatting');
    return
end

%make sure the intensities are positive
if sum(texttemp(:,4)<0)>0
    warndlg('Intensities must be positive values.','Incorrect Formatting');
    return
end

%make sure frame numbers are positive whole numbers
if or(sum(texttemp(:,5)<0)>0,sum(round(texttemp(:,5))~=texttemp(:,5))>0)
    warndlg('Frame number is a positive integer starting from 1.','Incorrect Formatting');
    return
end

if ncols==6
    if sum(texttemp(:,6)<0)>0
        warndlg('Labels are positive integers starting from 1.','Incorrect Formatting');
        return
    end
end

if and(ncols==5,get(handles.checkbox_stochastic,'Value'))
    set(handles.checkbox_stochastic,'Value',0)
    warndlg('To use stochastic bleaching, the input file must have a sixth "Label" column.','Stochastic Bleaching Invalid')
end

%check for sequential ordering
maxf=max(texttemp(:,5));

badorder=0;
for i=1:maxf
    if sum(texttemp(:,5)==i)<1
        badorder=1;
        break
    end
end

if badorder==1
    warndlg('Frames must be sequential whole numbers starting from 1.','Incorrect Formatting');
    return
end

%set the scale of the input data
if strcmp(button1,'nanometers')
    textput=texttemp;
else
    textput=texttemp;
    textput(:,1:3)=1000*textput(:,1:3);
end
close(h1);

%set output image size and frames from text file
sizestring=[num2str(size(textput,1)) ' points'];

Nfile=max(textput(:,5));
startim=1;
endim=Nfile;

basename=['output_' date '_' fname '.tif'];

set(handles.text_imfile,'String',num2str(Nfile))
set(handles.edit_startim,'String',num2str(startim))
set(handles.edit_endim,'String',num2str(endim))
set(handles.edit_fname,'String',basename)
set(handles.text_imsize,'String',sizestring)
set(handles.edit_slidenum,'String',num2str(startim));

if ncols==6
    set(handles.text_label,'String','File type:   labeled')
else
    set(handles.text_label,'String','File type:   unlabeled')
end

%display slider conditionally
if startim~=endim
    set(handles.slider1,'Enable','on')
    set(handles.slider1,'Min',startim);
    set(handles.slider1,'Max',endim);
    set(handles.slider1,'Value',startim);
    set(handles.text_slidemax,'Enable','on')
    set(handles.text_slidemin,'Enable','on')
    set(handles.text_slidemax,'String',num2str(endim));
    set(handles.text_slidemin,'String',num2str(startim));
else
    set(handles.slider1,'Enable','off')
    set(handles.text_slidemax,'Enable','off')
    set(handles.text_slidemin,'Enable','off')
end

pushbutton_center_Callback(hObject, eventdata, handles)
memupdate_Callback(hObject, eventdata, handles)
guidata(hObject,handles)

% --- Executes on button press in pushbutton_psf.
function pushbutton_psf_Callback(hObject, eventdata, handles)
global file2
global Npsf
global sx2
global sy2
global psfch

[fullname,cd_psf,~]=uigetfile({'*.tif','Image Stack (*.tif)'},'Pick a PSF input image stack.');

if fullname~=0
    prompt={'Enter the pixel calibration (nm/pixel): ','Enter the Z-spacing (nm/slice): '};
    name='PSF Calibration';
    numlines=1;
    defaultanswer={get(handles.edit_impxlsize,'String'),get(handles.edit_psfdz,'String')};
    answer=inputdlg(prompt,name,numlines,defaultanswer);
    
    if ~isempty(answer)
        set(handles.edit_impxlsize,'String',char(answer{1}))
        set(handles.edit_psfdz,'String',char(answer{2}))
        
        file2=[cd_psf fullname];
        
        fileinfo=imfinfo(file2);
        Npsf=length(fileinfo);
        
        [sy2, sx2]=size(imread(file2,1));
        sizestring=[num2str(sy2) ' x ' num2str(sx2)];
        
        if sy2~=sx2
            warndlg('The PSF images must be square in the X-Y plane.','PSF Size Invalid')
            return
        end
        
        set(handles.text_psfim,'String',num2str(Npsf))
        set(handles.text_psfsize,'String',sizestring)
        
        psfch=1;
        
        memupdate_Callback(hObject, eventdata, handles)
    end
end

% --- Executes on button press in pushbutton_iminfo.
function pushbutton_iminfo_Callback(hObject, eventdata, handles)
global Ikeep
global simch
global xpos
global ypos

fixq=get(handles.checkbox_fix,'Value');
startim=str2num(get(handles.edit_startim,'String'));

if simch==1
    sliderpos=str2num(get(handles.edit_slidenum,'String'));
    
    %get image
    if get(handles.checkbox_ram,'Value')
        Iout=Ikeep(:,:,sliderpos-startim+1);
    else
        try
            Iout=double(imread('BlurLab_temp.tif',sliderpos-startim+1));
        catch er1
            warndlg('The temporary stack file is missing -- please restart BlurLab and try your simulation again.','Missing Temp File')
            return
        end
    end
    
    axes(handles.axes1);
    if fixq==0
        imagesc(xpos/1000,ypos/1000,Iout)
    else
        intlow=str2num(get(handles.edit_lowfix,'String'));
        inthigh=str2num(get(handles.edit_highfix,'String'));
        imagesc(xpos/1000,ypos/1000,Iout,[intlow,inthigh])
    end
    axis equal
    set(gca,'color', [0.1 0.1 0.1])
    xlabel('x (um)')
    ylabel('y (um)')
    box on
    
    if get(handles.checkbox_histogram,'Value')
        histupdate_Callback(hObject, eventdata, handles, Iout)
    end
    
    axes(handles.axes1)
    impixelregion
else
    warndlg('No image information to examine.','No Images')
end

% --- Executes on button press in pushbutton_contrast.
function pushbutton_contrast_Callback(hObject, eventdata, handles)
global simch
if simch==1
    axes(handles.axes1)
    imcontrast(gca)
else
    warndlg('No image information to examine.','No Images')
end

% --- Executes on button press in pushbutton_play.
function pushbutton_play_Callback(hObject, eventdata, handles)
global simch
global Ikeep
global xpos
global ypos
global pvals
global intvals
global cd1
global basename
global cd_write

memq=get(handles.checkbox_ram,'Value');

if memq==1
    [sy,sx]=size(Ikeep(:,:,1));
else
    try
        [sy,sx]=size(imread('BlurLab_temp.tif',1));
    catch er1
        warndlg('The temporary disk file is missing -- please restart BlurLab and try your simulation again.','Sorry')
        return
    end
end

%get Z-stack parameters
zstackq=get(handles.checkbox_simZ,'Value');
minz=str2num(get(handles.edit_minz,'String'));
maxz=str2num(get(handles.edit_maxz,'String'));
dzstack=str2num(get(handles.edit_dzstack,'String'));
if zstackq==1
    zslice=linspace(minz,maxz,dzstack);
else
    zslice=str2num(get(handles.edit_zslice,'String'));
end

%handle cases with frapping and stochastic bleaching
if and(get(handles.checkbox_stochastic,'Value'),get(handles.checkbox_usenoise,'Value'))
    set(handles.checkbox_stochastic,'Value',0)
    warndlg('Stochastic bleaching cannot be used during playback.','Stochastic Bleaching Invalid')
    uiwait
end

if get(handles.checkbox_frap,'Value')
    set(handles.checkbox_frap,'Value',0)
    warndlg('FRAP cannot be used during playback.','FRAP Invalid')
    uiwait
end

if get(handles.checkbox_tirf,'Value')
    set(handles.checkbox_tirf,'Value',0)
    warndlg('TIRF cannot be used during playback.','TIRF Invalid')
    uiwait
end

writeq=get(handles.checkbox_write,'Value');
if writeq==1
    %make sure they want to save the playing file
    sureq = questdlg('Do you want to save this playback?', ...
        'Save Playback', ...
        'Yes', 'No', 'No');
    if strcmp(sureq,'No')
        writeq=0;
        set(handles.checkbox_write,'Value',0)
        set(handles.checkbox_savepars,'Value',0)
    else
        curr=get(handles.slider1,'Min');
        set(handles.slider1,'Value',curr)
        set(handles.edit_slidenum,'String',num2str(curr))
        set(handles.edit_znum,'String',num2str(1))
        set(handles.checkbox_loop,'Value',0)
        basename=get(handles.edit_fname,'String');
        if get(handles.checkbox_8bit,'Value')==1
            bitdepth=8;
        else
            bitdepth=16;
        end
    end
end

%assign cdout directory for image writing
if cd_write==0
    cdout=cd1;
else
    cdout=cd_write;
end

%check for redundant file output names and create unique, sequential names
if writeq==1
    l1=dir([cdout basename]);
    if ~isempty(l1)
        dotp=max(find(basename=='.'));
        numpos=max(findstr(basename,'_n'));
        if ~isempty(numpos)
            fnum=1+str2num(basename(numpos+2:dotp-1));
            basename=[basename(1:numpos-1) '_n' num2str(fnum) basename(dotp:end)];
        else
            basename=[basename(1:dotp-1) '_n1' basename(dotp:end)];
        end
        set(handles.edit_fname,'String',basename)
    end
end

%write output parameters file
if get(handles.checkbox_savepars,'Value')
    dotp=max(find(basename=='.'));
    parname=[basename(1:dotp) 'txt'];
    parametersave_Callback(hObject, eventdata, handles, cdout, parname)
end

%get boxcar parameters
boxcarq=get(handles.checkbox_boxcar,'Value');
loopstep=str2num(get(handles.edit_boxcar,'String'));

if simch==1
    if memq==1
        Nframe=size(Ikeep,3);
    else
        Nframe=length(imfinfo('BlurLab_temp.tif'));
    end
    
    set(handles.slider1,'Min',1)
    set(handles.slider1,'Max',Nframe)
    
    %play pause and reset movie
    startim=str2num(get(handles.edit_startim,'String'));
    endim=str2num(get(handles.edit_endim,'String'));
    
    if boxcarq==1
        endim=floor(endim/loopstep);
    end
    
    curr=str2num(get(handles.edit_slidenum,'String'));
    if curr==endim
        curr=startim;
    end
    
    lp=1;
    while lp==1
        q=curr-startim;       
        for i=curr:endim
            if and(get(handles.checkbox_stochastic,'Value'),get(handles.checkbox_usenoise,'Value'))
                set(handles.checkbox_stochastic,'Value',0)
            end
            
            q=q+1;
            if q>Nframe
                set(handles.edit_slidenum,'String',1);
                set(handles.text_slidemax,'String',num2str(Nframe));
                set(handles.text_slidemin,'String',num2str(1));
                set(handles.edit_endim,'String',num2str(Nframe));
                set(handles.edit_startim,'String',num2str(1));
                set(handles.slider1,'Min',1);
                set(handles.slider1,'Max',Nframe);
                set(handles.slider1,'Value',1);
                return
            else
                for p=1:length(zslice)
                    global pausech
                    if pausech==1
                        pausech=0;
                        return
                    end
                    
                    if memq==1
                        Iout=Ikeep(:,:,(q-1)*length(zslice)+p);
                    else
                        Iout=double(imread('BlurLab_temp.tif',(q-1)*length(zslice)+p));
                    end
                    
                    %turn off noise sources
                    if get(handles.checkbox_usenoise,'Value')
                        %scale intensity
                        Iout=Iout*str2num(get(handles.edit_sfc,'String'));
                        
                        %photobleaching
                        if get(handles.checkbox_bleach,'Value');
                            tau=str2num(get(handles.edit_tau,'String'));
                            Iout=Iout*exp(-(i-startim)/tau);
                        end
                        
                        %add Poisson/Shot noise
                        if get(handles.checkbox_poiss,'Value')
                            lam=str2num(get(handles.edit_poiss,'String'));
                            basal=str2num(get(handles.edit_basal,'String'));
                            Iout=poissimnoise(uint16(Iout),lam,basal);
                        end
                        
                        %add custom noise spectrum
                        if get(handles.checkbox_custom,'Value')
                            T1=custom_noise(pvals,1,outX*outY);
                            Iout=Iout+reshape(intvals(T1),outY,outX);
                        end
                        
                        %add Gaussian noise
                        if get(handles.checkbox_noise,'Value')
                            mu_gauss=str2num(get(handles.edit_mean,'String'));
                            std_gauss=str2num(get(handles.edit_std,'String'));
                            Iout=Iout+abs(normrnd(mu_gauss,std_gauss,sy,sx));
                        end
                    end
                    
                    axes(handles.axes1);
                    if get(handles.checkbox_fix,'Value')==0
                        imagesc(xpos/1000,ypos/1000,Iout)
                    else
                        intlow=str2num(get(handles.edit_lowfix,'String'));
                        inthigh=str2num(get(handles.edit_highfix,'String'));
                        imagesc(xpos/1000,ypos/1000,Iout,[intlow,inthigh])
                    end
                    axis equal
                    set(gca,'color', [0.1 0.1 0.1])
                    xlabel('x (um)')
                    ylabel('y (um)')
                    box on
                    
                    if get(handles.checkbox_histogram,'Value')
                        histupdate_Callback(hObject, eventdata, handles, Iout)
                    end
                    
                    %write output images
                    if writeq==1
                        er1=0;
                        while er1==0
                            try
                                if bitdepth==8
                                    imwrite(uint8(Iout),[cdout basename],'Compression','none','Writemode','append');
                                else
                                    imwrite(uint16(Iout),[cdout basename],'Compression','none','Writemode','append');
                                end
                                er1=1;
                            catch
                                pause(rand)
                            end
                        end
                    end
                    
                    set(handles.edit_znum,'String',num2str(p));
                    set(handles.edit_slidenum,'String',num2str(i));
                    set(handles.slider1,'Value',i);
                    set(handles.text_intmin,'String',num2str(1/100*round(100*min(min(Iout)))));
                    set(handles.text_intmax,'String',num2str(1/100*round(100*max(max(Iout)))));
                    pause(1/str2num(get(handles.edit_fps,'String')))
                end
            end
        end
        %handling playback looping
        if get(handles.checkbox_loop,'Value')
            lp=1;
            curr=startim;
        else
            lp=0;
        end
    end
else
    warndlg('There is no data to play back.','Nothing to play back')
end

% --- Executes on button press in pushbutton_loadstk.
function pushbutton_loadstk_Callback(hObject, eventdata, handles)
global Ikeep
global xpos
global ypos
global loadch
global psfch
global simch
global ncols

[fullname,cd1,~]=uigetfile({'*.tif','Image Stack (*.tif)'},'Choose an image stack for playback.');

memupdate_Callback(hObject, eventdata, handles)

if fullname~=0
    %handle cases with frapping and stochastic bleaching
    if and(get(handles.checkbox_bleach,'Value'),get(handles.checkbox_usenoise,'Value'))
        set(handles.checkbox_stochastic,'Value',0)
    end
    
    file3=[cd1 fullname];
    
    %get image data
    fileinfo=imfinfo(file3);
    Nfile=length(fileinfo);
    I0=imread(file3,1);
    
    %turn off z-stacking
    set(handles.checkbox_simZ,'Value',0);
    set(handles.edit_znum,'String',num2str(1));
    
    %check channels
    if length(size(I0))>2
        warndlg('Each image in the stack should have a single channel.','File Format Error')
        return
    end
    
    loadch=1;
    psfch=1;
    simch=1;
    ncols=5;
    
    %change basename
    set(handles.edit_fname,'String',['output_' date '_'  fullname]);
    
    startim=1;
    endim=Nfile;
    
    %get image and stack size
    [sy, sx]=size(imread(file3,1));
    sizestring=[num2str(sy) ' x ' num2str(sx)];
    set(handles.edit_outX,'String',num2str(sx));
    set(handles.edit_outY,'String',num2str(sy));
    set(handles.text_imsize,'String',sizestring);
    set(handles.text_imfile,'String',num2str(Nfile));
    
    %calibrate image
    imcalb=str2num(get(handles.edit_impxlsize,'String'));
    xpos=0:imcalb:(sx-1)*imcalb;
    ypos=0:imcalb:(sy-1)*imcalb;
    
    %setup slider
    set(handles.slider1,'Min',startim);
    set(handles.slider1,'Max',endim);
    set(handles.slider1,'Value',startim);
    set(handles.text_slidemax,'String',num2str(endim));
    set(handles.text_slidemin,'String',num2str(startim));
    set(handles.edit_slidenum,'String',num2str(startim));
    set(handles.edit_endim,'String',num2str(endim));
    set(handles.edit_startim,'String',num2str(startim));
    
    if get(handles.checkbox_ram,'Value')
        Ikeep=zeros(size(I0,1),size(I0,2),Nfile);
        for i=1:Nfile
            Ikeep(:,:,i)=imread(file3,i);
        end
    else
        if ~isempty(dir('BlurLab_temp.tif'))
            delete BlurLab_temp.tif
        end
        copyfile(file3,'BlurLab_temp.tif')
    end
    
    memupdate_Callback(hObject, eventdata, handles)
    guidata(hObject,handles)
end

% --- Executes on button press in pushbutton_boxselect.
function pushbutton_boxselect_Callback(hObject, eventdata, handles)
global psfch
global loadch
global frapch

if and(psfch==1,loadch==1)
    axes(handles.axes1);
    rect1=getrect(gca);
    
    ans1=0;
    while ans1==0
        hold on
        p1=plot([rect1(1),rect1(1)+rect1(3)],[rect1(2),rect1(2)],'w--','Linewidth',2);
        p2=plot([rect1(1),rect1(1)+rect1(3)],[rect1(2)+rect1(4),rect1(2)+rect1(4)],'w--','Linewidth',2);
        p3=plot([rect1(1),rect1(1)],[rect1(2),rect1(2)+rect1(4)],'w--','Linewidth',2);
        p4=plot([rect1(1)+rect1(3),rect1(1)+rect1(3)],[rect1(2),rect1(2)+rect1(4)],'w--','Linewidth',2);
        hold off
        
        ans2=questdlg('Are you satisfied with the selected region?','FRAP Region Accept',...
            'Yes','No','Cancel','Yes');
        if strcmp('Yes',ans2)
            ans1=1;
            global xrectmin
            global xrectmax
            global yrectmin
            global yrectmax
            
            xrectmin=1000*(rect1(1));
            xrectmax=1000*(rect1(1)+rect1(3));
            yrectmin=1000*(rect1(2));
            yrectmax=1000*(rect1(2)+rect1(4));
            delete(p1)
            delete(p2)
            delete(p3)
            delete(p4)
            frapch=1;
            set(handles.checkbox_frap,'Value',1)
            return
        elseif strcmp('No',ans2)
            delete(p1)
            delete(p2)
            delete(p3)
            delete(p4)
            rect1=getrect(gca);
        else
            delete(p1)
            delete(p2)
            delete(p3)
            delete(p4)
            ans1=1;
            set(handles.checkbox_frap,'Value',0)
            return
        end
    end
else
    warndlg('To select a FRAP region, first load input data and a PSF file.','Load Files for FRAP')
    return
end

% --- Executes on button press in pushbutton_custom.
function pushbutton_custom_Callback(hObject, eventdata, handles)
global pvals
global intvals
global noisech

[fullname2,cd2,filter2]=uigetfile({'*.tif','Image Stack (*.tif)';'*.txt','Text Input (*.txt)'},'Pick a custom noise file.');

memupdate_Callback(hObject, eventdata, handles)

if fullname2~=0
    [fn]=max(find(fullname2=='.'));
    fname2=fullname2(1:fn-1);
    file2=[cd2 fullname2];
    
    %convert tif to text file
    if filter2==1
        fileinfo2=imfinfo(file2);
        Nnoise=length(fileinfo2);
        
        %create output text file for image spots
        txtname2=['text_noise_' fname2 '_' date '.txt'];
        fid2 = fopen([cd2 txtname2],'w');
        
        %*****TEXT OUT******
        %establish max and min noise values
        minnoise=2^16-1;
        maxnoise=0;
        for i=1:Nnoise
            Inoise=double(imread(file2,i));
            maxtemp=max(max(Inoise));
            mintemp=min(min(Inoise));
            
            if maxtemp>maxnoise
                maxnoise=maxtemp;
            end
            if mintemp<minnoise
                minnoise=mintemp;
            end
        end
        
        %establish intensity vales
        intvals=minnoise:maxnoise;
        
        %construct stack histogram
        pvals=zeros(1,length(intvals));
        for i=1:Nnoise
            Inoise=imread(file2,i);
            
            [pvals_temp,~]=hist(Inoise(Inoise>=0),intvals);
            
            pvals=pvals+pvals_temp;
        end
        pvals=pvals/sum(pvals);
        
        fprintf(fid2,'%6f %6.4f\n',[intvals;pvals]);
        fclose(fid2);
        noisech=1;
    else
        textnoise=load(file2);
        
        if size(textnoise,2)~=2
            warndlg('This file should have two columns.','Incorrect Format')
            return
        end
        
        if sum(sum(textnoise<0))>0
            warndlg('All values in this file should be positive.','Incorrect Format')
            return
        end
        
        intvals=textnoise(:,1);
        pvals=textnoise(:,2)/sum(textnoise(:,2));
        noisech=1;
    end
end
memupdate_Callback(hObject, eventdata, handles)
guidata(hObject,handles)


% --- Executes on button press in pushbutton_center.
function pushbutton_center_Callback(hObject, eventdata, handles)
global textput
global sy
global sx
global psfch
global loadch

if loadch==1
pxcalb=str2num(get(handles.edit_impxlsize,'String'));
bord=str2num(get(handles.edit_border,'String'));

xmin=min(textput(:,1));
xmax=max(textput(:,1));
xmid=(xmin+xmax)/2;
ymin=min(textput(:,2));
ymax=max(textput(:,2));
ymid=(ymin+ymax)/2;

sx=round((2*bord+xmax-xmin)/pxcalb);
sy=round((2*bord+ymax-ymin)/pxcalb);

xmidsx=sx*pxcalb/2;
ymidsy=sy*pxcalb/2;

offx=xmidsx-xmid;
offy=ymidsy-ymid;

set(handles.edit_outX,'String',num2str(sx))
set(handles.edit_outY,'String',num2str(sy))
set(handles.edit_offsetx,'String',num2str(offx))
set(handles.edit_offsety,'String',num2str(offy))
else
    warndlg('No input data to center.','Load Input File')
end
% if get(handles.checkbox_flip,'Value')
%     set(handles.edit_outX,'String',num2str(sy))
%     set(handles.edit_outY,'String',num2str(sx))
%     set(handles.edit_offsetx,'String',num2str(offy))
%     set(handles.edit_offsety,'String',num2str(offx))
% else
%     set(handles.edit_outX,'String',num2str(sx))
%     set(handles.edit_outY,'String',num2str(sy))
%     set(handles.edit_offsetx,'String',num2str(offx))
%     set(handles.edit_offsety,'String',num2str(offy))
% end


% --- Executes on button press in pushbutton_directory.
function pushbutton_directory_Callback(hObject, eventdata, handles)
global cd_write
[cd_write]=uigetdir;
if cd_write==0
else
    cd_write=[cd_write '\'];
end


% --- Executes on button press in pushbutton_generate.
function pushbutton_generate_Callback(hObject, eventdata, handles)
%get parameters from gui
nframes=str2num(get(handles.edit_randframes,'String'));
meanI=str2num(get(handles.edit_meanint,'String'));
D=str2num(get(handles.edit_movestd,'String'));
npts=str2num(get(handles.edit_points,'String'));
Lx=str2num(get(handles.edit_Lx,'String'));
Ly=str2num(get(handles.edit_Ly,'String'));
Lz=str2num(get(handles.edit_Lz,'String'));

random_box(nframes,npts,meanI,D,Lx,Ly,Lz)
%pushbutton_load_Callback(hObject, eventdata, handles)


% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
%helppopup('BlurLab_Help.txt','BlurLab Quick Guide')
try
    fname=['.\documentation\BlurLab_v0.9_User_Manual.pdf'];
    open(fname)
catch
    try
        fname=['./documentation/BlurLab_v0.9_User_Manual.pdf'];
        open(fname)
    catch
        warndlg('Could not locate help file in Documentation folder.','Missing Help File')
    end
end

% --- Executes on button press in pushbutton_stop.
function pushbutton_stop_Callback(hObject, eventdata, handles)
global stopch
stopch=1;

% --- Executes on button press in pushbutton_pause.
function pushbutton_pause_Callback(hObject, eventdata, handles)
global pausech
pausech=1;

% --- Executes on button press in pushbutton_clear.
function pushbutton_clear_Callback(hObject, eventdata, handles)
clearvars
close all
%clear temp image file
if ~isempty(dir('BlurLab_temp.tif'))
    delete BlurLab_temp.tif
end
pause(0.1)
BlurLab_3D

% --- Executes on button press in pushbutton_exit.
function pushbutton_exit_Callback(hObject, eventdata, handles)
%clear temp image file
if ~isempty(dir('BlurLab_temp.tif'))
    delete BlurLab_temp.tif
end
clearvars
close all

%****************SLIDER*************************************************
%***********************************************************************
% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global loadch
global psfch
global simch
global Ikeep
global xpos
global ypos

fixq=get(handles.checkbox_fix,'Value');
boxcarq=get(handles.checkbox_boxcar,'Value');
loopstep=str2num(get(handles.edit_boxcar,'String'));

if and(loadch==1,psfch==1)
    startim=str2num(get(handles.edit_startim,'String'));
    endim=str2num(get(handles.edit_endim,'String'));
    
    if or(and(startim==endim,boxcarq==0)...
            ,and(boxcarq==1,startim==floor(endim/loopstep)))
        set(handles.slider1,'Enable','off')
        set(handles.text_slidemax,'Enable','off')
        set(handles.text_slidemin,'Enable','off')
        set(handles.slider1,'Min',startim)
        set(handles.slider1,'Max',startim)
        set(handles.slider1,'Value',startim)
    else
        set(handles.slider1,'Enable','on')
        set(handles.text_slidemax,'Enable','on')
        set(handles.text_slidemin,'Enable','on')
        
        if boxcarq==0
            set(handles.slider1,'Min',startim)
            set(handles.slider1,'Max',endim)
            set(handles.text_slidemax,'String',num2str(endim));
            set(handles.text_slidemin,'String',num2str(startim));
        else
            set(handles.slider1,'Min',startim)
            set(handles.slider1,'Max',floor(endim/loopstep))
            set(handles.text_slidemax,'String',num2str(floor(endim/loopstep)));
            set(handles.text_slidemin,'String',num2str(startim));
        end
        
        %sliderpos=round(str2num(get(handles.edit_slidenum,'String')));
        sliderpos=round(get(handles.slider1,'Value'));
        set(handles.edit_slidenum,'String',num2str(sliderpos));
        
        if simch==1
            if get(handles.checkbox_simZ,'Value')==0
                q=sliderpos-startim+1;
            else
                fnum=str2num(get(handles.edit_slidenum,'String'));
                minz=str2num(get(handles.edit_minz,'String'));
                maxz=str2num(get(handles.edit_maxz,'String'));
                dzstack=str2num(get(handles.edit_dzstack,'String'));
                
                znum=str2num(get(handles.edit_znum,'String'));
                znum_max=floor((maxz-minz)/dzstack)+1;
                
                q=(fnum-startim)*znum_max+znum;
            end
            
            if q<1
                q=1;
            end
            
            if get(handles.checkbox_ram,'Value')
                cutq=size(Ikeep,3);
            else
                cutq=length(imfinfo('BlurLab_temp.tif'));
            end
            
            if cutq>=q
                if get(handles.checkbox_ram,'Value')
                    Iout=Ikeep(:,:,q);
                else
                    Iout=double(imread('BlurLab_temp.tif',q));
                end
                
                %make output plot
                axes(handles.axes1);
                if fixq==0
                    imagesc(xpos/1000,ypos/1000,Iout)
                else
                    intlow=str2num(get(handles.edit_lowfix,'String'));
                    inthigh=str2num(get(handles.edit_highfix,'String'));
                    imagesc(xpos/1000,ypos/1000,Iout,[intlow,inthigh])
                end
                axis equal
                set(gca,'color', [0.1 0.1 0.1])
                xlabel('x (um)')
                ylabel('y (um)')
                box on
                
                if get(handles.checkbox_histogram,'Value')
                    histupdate_Callback(hObject, eventdata, handles, Iout)
                end
            else
                set(handles.slider1,'Value',get(handles.slider1,'Min'))
                set(handles.edit_slidenum,'String',num2str(get(handles.slider1,'Min')))
                warndlg('Slider value out of bounds of current image stack.','Slider Out of Bounds')
                return
            end
        end
    end
    guidata(hObject,handles)
else
    warndlg('Load input and PSF files to activate.','Slider Inactive')
end
% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

%******************EDIT BOXES*******************************************
%***********************************************************************
function edit_mean_Callback(hObject, eventdata, handles)
global mu_gauss
mu_gauss=str2num(get(handles.edit_mean,'String'));
if mu_gauss<1
    mu_gauss=1;
    set(handles.edit_mean,'String',num2str(mu_gauss))
end
% --- Executes during object creation, after setting all properties.
function edit_mean_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_startim_Callback(hObject, eventdata, handles)
startim=round(str2num(get(handles.edit_startim,'String')));
endim=round(str2num(get(handles.edit_endim,'String')));

if startim<1
    startim=1;
    set(handles.edit_startim,'String','1')
elseif startim>endim
    endim=startim;
    set(handles.edit_endim,'String',num2str(endim))
end

slider1_Callback(hObject, eventdata, handles)
set(handles.slider1,'Value',startim)
set(handles.edit_startim,'String',num2str(startim))
set(handles.edit_endim,'String',num2str(endim))
set(handles.edit_slidenum,'String',num2str(startim))
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function edit_startim_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_endim_Callback(hObject, eventdata, handles)
global Nfile
startim=round(str2num(get(handles.edit_startim,'String')));
endim=round(str2num(get(handles.edit_endim,'String')));

if endim>Nfile
    endim=Nfile;
    set(handles.edit_endim,'String',num2str(Nfile))
elseif or(endim<startim,endim<1)
    endim=startim;
    set(handles.edit_endim,'String',num2str(endim))
end

slider1_Callback(hObject, eventdata, handles)
set(handles.slider1,'Value',startim)
set(handles.edit_startim,'String',num2str(startim))
set(handles.edit_endim,'String',num2str(endim))
set(handles.edit_slidenum,'String',num2str(startim))
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function edit_endim_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_fname_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function edit_fname_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_slidenum_Callback(hObject, eventdata, handles)
global simch
global Ikeep
global xpos
global ypos

fixq=get(handles.checkbox_fix,'Value');
startim=str2num(get(handles.edit_startim,'String'));

sliderpos=round(str2num(get(handles.edit_slidenum,'String')));
set(handles.edit_slidenum,'String',num2str(sliderpos))

if sliderpos>get(handles.slider1,'Max')
    set(handles.slider1,'Value',get(handles.slider1,'Max'))
    set(handles.edit_slidenum,'String',get(handles.slider1,'Max'))
elseif sliderpos<get(handles.slider1,'Min')
    set(handles.slider1,'Value',get(handles.slider1,'Min'))
    set(handles.edit_slidenum,'String',get(handles.slider1,'Min'))
else
    set(handles.slider1,'Value',sliderpos)
end

sliderpos=round(str2num(get(handles.edit_slidenum,'String')));

if simch==1
    if get(handles.checkbox_simZ,'Value')==0
        q=sliderpos-startim+1;
    else
        fnum=str2num(get(handles.edit_slidenum,'String'));
        minz=str2num(get(handles.edit_minz,'String'));
        maxz=str2num(get(handles.edit_maxz,'String'));
        dzstack=str2num(get(handles.edit_dzstack,'String'));
        
        znum=str2num(get(handles.edit_znum,'String'));
        znum_max=floor((maxz-minz)/dzstack)+1;
        
        q=(fnum-startim)*znum_max+znum;
    end
    
    if get(handles.checkbox_ram,'Value')
        Iout=Ikeep(:,:,q);
    else
        Iout=double(imread('BlurLab_temp.tif',q));
    end
    
    %get image size
    [sy,sx]=size(Iout);
    
    %turn off noise sources
    if get(handles.checkbox_usenoise,'Value')==1
        %scale intensity
        Iout=Iout*str2num(get(handles.edit_sfc,'String'));
        
        %photobleaching
        if get(handles.checkbox_bleach,'Value');
            tau=str2num(get(handles.edit_tau,'String'));
            Iout=Iout*exp(-(i-startim)/tau);
        end
        
        %add Poisson/Shot noise
        if get(handles.checkbox_poiss,'Value')
            lam=str2num(get(handles.edit_poiss,'String'));
            basal=str2num(get(handles.edit_basal,'String'));
            Iout=poissimnoise(uint16(Iout),lam,basal);
        end
        
        %add custom noise spectrum
        if get(handles.checkbox_custom,'Value')
            T1=custom_noise(pvals,1,outX*outY);
            Iout=Iout+reshape(intvals(T1),outY,outX);
        end
        
        %add Gaussian noise
        if get(handles.checkbox_noise,'Value')
            mu_gauss=str2num(get(handles.edit_mean,'String'));
            std_gauss=str2num(get(handles.edit_std,'String'));
            Iout=Iout+abs(normrnd(mu_gauss,std_gauss,sy,sx));
        end
    end
    
    %make output plot
    axes(handles.axes1);
    if fixq==0
        imagesc(xpos/1000,ypos/1000,Iout)
    else
        intlow=str2num(get(handles.edit_lowfix,'String'));
        inthigh=str2num(get(handles.edit_highfix,'String'));
        imagesc(xpos/1000,ypos/1000,Iout,[intlow,inthigh])
    end
    axis equal
    set(gca,'color', [0.1 0.1 0.1])
    xlabel('x (um)')
    ylabel('y (um)')
    box on
    
    if get(handles.checkbox_histogram,'Value')
        histupdate_Callback(hObject, eventdata, handles, Iout)
    end
end
guidata(hObject,handles)
pause(0.01)
% --- Executes during object creation, after setting all properties.
function edit_slidenum_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_psfdz_Callback(hObject, eventdata, handles)
global dz
dz=str2num(get(handles.edit_psfdz,'String'));
if dz<=1
    dz=10;
    set(handles.edit_psfdz,'String',num2str(dz))
end
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function edit_psfdz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_impxlsize_Callback(hObject, eventdata, handles)
global imcalb
global loadch
imcalb=str2num(get(handles.edit_impxlsize,'String'));
if imcalb<=0
    imcalb=80;
    set(handles.edit_impxlsize,'String',num2str(imcalb))
end

if loadch==1
    pushbutton_center_Callback(hObject, eventdata, handles)
end
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function edit_impxlsize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_zslice_Callback(hObject, eventdata, handles)
global psfch
global Npsf
global dz
if psfch==0
    warndlg('Load PSF file, then select Z slice.','Load PSF')
    set(handles.edit_zslice,'String' ,'0')
end
% --- Executes during object creation, after setting all properties.
function edit_zslice_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_outX_Callback(hObject, eventdata, handles)
global sx
global loadch

if loadch==0
    set(handles.edit_outX,'String','0');
    warndlg('Load input file, then set output size.','Load File')
else
    outX=ceil(str2num(get(handles.edit_outX,'String')));
    if outX<=1
        outX=sx;
        set(handles.edit_outX,'String',num2str(outX));
    end
    guidata(hObject,handles)
end
% --- Executes during object creation, after setting all properties.
function edit_outX_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_outY_Callback(hObject, eventdata, handles)
global sy
global loadch

if loadch==0
    set(handles.edit_outY,'String','0');
    warndlg('Load input file, then set output size.','Load File')
else
    outY=ceil(str2num(get(handles.edit_outY,'String')));
    if outY<=1
        outY=sy;
        set(handles.edit_outY,'String',num2str(outY));
    end
    guidata(hObject,handles)
end
% --- Executes during object creation, after setting all properties.
function edit_outY_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_points_Callback(hObject, eventdata, handles)
npoints=ceil(str2num(get(handles.edit_points,'String')));
if npoints<1
    npoints=1;
    set(handles.edit_points,'String',num2str(npoints))
end
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function edit_points_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_meanint_Callback(hObject, eventdata, handles)
global meanint
meanint=ceil(str2num(get(handles.edit_meanint,'String')));
if meanint<1
    meanint=1;
    set(handles.edit_meanint,'String',num2str(meanint))
end
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function edit_meanint_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_movestd_Callback(hObject, eventdata, handles)
global movestd
movestd=ceil(str2num(get(handles.edit_movestd,'String')));
if movestd<1
    movestd=2;
    set(handles.edit_movestd,'String',num2str(movestd))
end
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function edit_movestd_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_sfc_Callback(hObject, eventdata, handles)
sfc=str2num(get(handles.edit_sfc,'String'));
if sfc<0
    sfc=1;
    set(handles.edit_sfc,'String',num2str(sfc))
end
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function edit_sfc_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_std_Callback(hObject, eventdata, handles)
global std_gauss
std_gauss=str2num(get(handles.edit_std,'String'));
if std_gauss<1
    std_gauss=1;
    set(handles.edit_std,'String',num2str(std_gauss))
end
% --- Executes during object creation, after setting all properties.
function edit_std_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_poiss_Callback(hObject, eventdata, handles)
if or(str2num(get(handles.edit_poiss,'String'))<0,isempty(str2num(get(handles.edit_poiss,'String'))))
    set(handles.edit_poiss,'String',num2str(1))
end
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function edit_poiss_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_lowfix_Callback(hObject, eventdata, handles)
intlow=str2num(get(handles.edit_lowfix,'String'));
inthigh=str2num(get(handles.edit_highfix,'String'));
if intlow<1
    set(handles.edit_lowfix,'String',num2str(1))
end
if or(inthigh<intlow,isempty(intlow))
    set(handles.edit_highfix,'String',num2str(intlow+1))
end
edit_slidenum_Callback(hObject, eventdata, handles)
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function edit_lowfix_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_highfix_Callback(hObject, eventdata, handles)
intlow=str2num(get(handles.edit_lowfix,'String'));
inthigh=str2num(get(handles.edit_highfix,'String'));
if inthigh>(2^16-1)
    set(handles.edit_highfix,'String',num2str(2^16-1))
end
if or(inthigh<intlow,isempty(inthigh))
    set(handles.edit_lowfix,'String',num2str(inthigh-1))
end
edit_slidenum_Callback(hObject, eventdata, handles)
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function edit_highfix_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_minz_Callback(hObject, eventdata, handles)
global Npsf
global psfch
minz=str2num(get(handles.edit_minz,'String'));
maxz=str2num(get(handles.edit_maxz,'String'));
dz=str2num(get(handles.edit_psfdz,'String'));

if psfch==0
    warndlg('Load PSF file before setting Z bounds.','Load PSF.')
    return
end

if minz>(maxz-1)
    set(handles.edit_minz,'String',num2str(maxz-1))
end
edit_dzstack_Callback(hObject, eventdata, handles)
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function edit_minz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_maxz_Callback(hObject, eventdata, handles)
global Npsf
global psfch
minz=str2num(get(handles.edit_minz,'String'));
maxz=str2num(get(handles.edit_maxz,'String'));
dz=str2num(get(handles.edit_psfdz,'String'));

if psfch==0
    warndlg('Load PSF file before setting Z bounds.','Load PSF.')
    return
end

if maxz<(minz+1)
    set(handles.edit_maxz,'String',num2str(minz+1))
end
edit_dzstack_Callback(hObject, eventdata, handles)
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function edit_maxz_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_dzstack_Callback(hObject, eventdata, handles)
minZ=str2num(get(handles.edit_minz,'String'));
maxZ=str2num(get(handles.edit_maxz,'String'));
dzstack=str2num(get(handles.edit_dzstack,'String'));
psfdz=str2num(get(handles.edit_psfdz,'String'));

if dzstack>floor((maxZ-minZ)/psfdz)
    set(handles.edit_dzstack,'String',num2str(floor((maxZ-minZ)/psfdz)));
end

if or(dzstack<1,length(dzstack)==0)
    set(handles.edit_dzstack,'String','1');
end
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function edit_dzstack_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_fps_Callback(hObject, eventdata, handles)
fps=str2num(get(handles.edit_fps,'String'));
if or(or(fps>100,fps<0.1),isempty(fps))
    set(handles.edit_fps,'String',num2str(1))
end
guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function edit_fps_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_frapframe_Callback(hObject, eventdata, handles)
startim=str2num(get(handles.edit_startim,'String'));
endim=str2num(get(handles.edit_endim,'String'));
frapframe=str2num(get(handles.edit_frapframe,'String'));

if or(frapframe<startim,frapframe<1)
    frapframe=startim;
elseif frapframe>endim
    frapframe=endim;
end
set(handles.edit_frapframe,'String',num2str(frapframe));
% --- Executes during object creation, after setting all properties.
function edit_frapframe_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_basal_Callback(hObject, eventdata, handles)
if str2num(get(handles.edit_basal,'String'))<0
    set(handles.edit_basal,'String',num2str(0))
end
% --- Executes during object creation, after setting all properties.
function edit_basal_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_frap_Callback(hObject, eventdata, handles)
frapper=str2num(get(handles.edit_frap,'String'));
if frapper<0
    set(handles.edit_frap,'String','0')
    frapper=0;
elseif frapper>100
    set(handles.edit_frap,'String','100')
    frapper=100;
end
if length(frapper)==0
    set(handles.edit_frap,'String','100')
    frapper=100;
end
% --- Executes during object creation, after setting all properties.
function edit_frap_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_znum_Callback(hObject, eventdata, handles)
global Ikeep
global simch
global xpos
global ypos
global pvals
global intvals

sfc=str2num(get(handles.edit_sfc,'String'));
startim=str2num(get(handles.edit_startim,'String'));

fixq=get(handles.checkbox_fix,'Value');
zstackq=get(handles.checkbox_simZ,'Value');

outX=str2num(get(handles.edit_outX,'String'));
outY=str2num(get(handles.edit_outY,'String'));

if zstackq==0
    set(handles.edit_znum,'String',num2str(1));
else
    if simch==1
        fnum=str2num(get(handles.edit_slidenum,'String'));
        minz=str2num(get(handles.edit_minz,'String'));
        maxz=str2num(get(handles.edit_maxz,'String'));
        dzstack=str2num(get(handles.edit_dzstack,'String'));
        
        znum=str2num(get(handles.edit_znum,'String'));
        znum_max=floor((maxz-minz)/dzstack)+1;
        
        %handle errors in z slice selection
        if znum<1
            znum=1;
            set(handles.edit_znum,'String',num2str(znum));
        elseif znum>znum_max
            znum=znum_max;
            set(handles.edit_znum,'String',num2str(znum));
        end
        
        q=(fnum-startim)*znum_max+znum;
        
        if get(handles.checkbox_ram,'Value')
            Iout=Ikeep(:,:,q);
        else
            Iout=double(imread('BlurLab_temp.tif',q));
        end
        
        if get(handles.checkbox_usenoise,'Value')
            %scale intensity
            Iout=Iout*str2num(get(handles.edit_sfc,'String'));
            
            %add custom noise spectrum
            if get(handles.checkbox_custom,'Value')
                T1=custom_noise(pvals,1,outX*outY);
                Iout=Iout+reshape(intvals(T1),outY,outX);
            end
            
            %add Poisson/Shot noise
            if get(handles.checkbox_poiss,'Value')
                lam=str2num(get(handles.edit_poiss,'String'));
                basal=str2num(get(handles.edit_basal,'String'));
                Iout=poissimnoise(uint16(Iout),lam,basal);
            end
            
            if get(handles.checkbox_noise,'Value')
                mu_gauss=str2num(get(handles.edit_mean,'String'));
                std_gauss=str2num(get(handles.edit_std,'String'));
                Iout=Iout+abs(normrnd(mu_gauss,std_gauss,outY,outX));
            end
        end
        
        %make output plot
        axes(handles.axes1);
        if fixq==0
            imagesc(xpos/1000,ypos/1000,Iout)
        else
            intlow=str2num(get(handles.edit_lowfix,'String'));
            inthigh=str2num(get(handles.edit_highfix,'String'));
            imagesc(xpos/1000,ypos/1000,Iout,[intlow,inthigh])
        end
        axis equal
        set(gca,'color', [0.1 0.1 0.1])
        xlabel('x (um)')
        ylabel('y (um)')
        box on
        
        if get(handles.checkbox_histogram,'Value')
            histupdate_Callback(hObject, eventdata, handles, Iout)
        end
    end
    guidata(hObject,handles)
    pause(0.01)
end
% --- Executes during object creation, after setting all properties.
function edit_znum_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_boxcar_Callback(hObject, eventdata, handles)
boxcar=str2num(get(handles.edit_boxcar,'String'));
startim=str2num(get(handles.edit_startim,'String'));
endim=str2num(get(handles.edit_endim,'String'));

if boxcar<1
    set(handles.edit_boxcar,'String','1')
elseif boxcar>(endim-startim+1)
    set(handles.edit_boxcar,'String',num2str(endim-startim+1))
end
edit_endim_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function edit_boxcar_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_tau_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function edit_tau_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_offsetx_Callback(hObject, eventdata, handles)
ofx=str2num(get(handles.edit_offsetx,'String'));
set(handles.edit_offsetx,'String',num2str(round(ofx)))
% --- Executes during object creation, after setting all properties.
function edit_offsetx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_offsety_Callback(hObject, eventdata, handles)
ofy=str2num(get(handles.edit_offsety,'String'));
set(handles.edit_offsety,'String',num2str(round(ofy)))
% --- Executes during object creation, after setting all properties.
function edit_offsety_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_framecalb_Callback(hObject, eventdata, handles)
% --- Executes during object creation, after setting all properties.
function edit_framecalb_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_Lx_Callback(hObject, eventdata, handles)
if str2num(get(handles.edit_Lx,'String'))<0.05
    set(handles.edit_Lx,'String',num2str(0.05))
end
% --- Executes during object creation, after setting all properties.
function edit_Lx_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_Ly_Callback(hObject, eventdata, handles)
if str2num(get(handles.edit_Ly,'String'))<0.05
    set(handles.edit_Ly,'String',num2str(0.05))
end
% --- Executes during object creation, after setting all properties.
function edit_Ly_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_Lz_Callback(hObject, eventdata, handles)
if str2num(get(handles.edit_Lz,'String'))<0.05
    set(handles.edit_Lz,'String',num2str(0.05))
end
% --- Executes during object creation, after setting all properties.
function edit_Lz_CreateFcn(hObject, eventdata, handles)


function edit_randframes_Callback(hObject, eventdata, handles)
if str2num(get(handles.edit_randframes,'String'))<1
    set(handles.edit_randframes,'String',num2str(1))
end
% --- Executes during object creation, after setting all properties.
function edit_randframes_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_tirfdecay_Callback(hObject, eventdata, handles)
if str2num(get(handles.edit_tirfdecay,'String'))<1
    set(handles.edit_tirfdecay,'String',num2str(1))
end 
% --- Executes during object creation, after setting all properties.
function edit_tirfdecay_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_tirfZ_Callback(hObject, eventdata, handles)
global loadch
if loadch==1
    global textput
    if str2num(get(handles.edit_tirfZ,'String'))>min(textput(:,3))
        set(handles.edit_tirfZ,'String',num2str(min(textput(:,3))))
    end
end
% --- Executes during object creation, after setting all properties.
function edit_tirfZ_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_border_Callback(hObject, eventdata, handles)
pxcalb=str2num(get(handles.edit_impxlsize,'String'));
if str2num(get(handles.edit_border,'String'))<3*pxcalb
    set(handles.edit_border,'String',num2str(3*pxcalb))
end
% --- Executes during object creation, after setting all properties.
function edit_border_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%**************CHECKBOXES*************************************
%*************************************************************

% --- Executes on button press in checkbox_usenoise.
function checkbox_usenoise_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_boxcar.
function checkbox_boxcar_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_tirf.
function checkbox_tirf_Callback(hObject, eventdata, handles)
global loadch
if loadch==1
    global textput
    set(handles.edit_tirfZ,'String',num2str(min(textput(:,3))))
end

% --- Executes on button press in checkbox_write.
function checkbox_write_Callback(hObject, eventdata, handles)
if get(handles.checkbox_write,'Value')
    set(handles.checkbox_savepars,'Value',1)
else
    set(handles.checkbox_savepars,'Value',0)
end

% --- Executes on button press in checkbox_psfnorm.
function checkbox_psfnorm_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_bleach.
function checkbox_bleach_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_flip.
function checkbox_flip_Callback(hObject, eventdata, handles)
global loadch
if loadch==1
    pushbutton_center_Callback(hObject, eventdata, handles)
end

% --- Executes on button press in checkbox_loop.
function checkbox_loop_Callback(hObject, eventdata, handles)

% --- Executes on button press in checkbox_stochastic.
function checkbox_stochastic_Callback(hObject, eventdata, handles)
global ncols
if ncols~=6
    set(handles.checkbox_stochastic,'Value',0)
    warndlg('To use stochastic bleaching, the input file must have a sixth "Label" column.','Stochastic Bleaching Invalid')
end

% --- Executes on button press in checkbox_savepars.
function checkbox_savepars_Callback(hObject, eventdata, handles)
if get(handles.checkbox_savepars,'Value')
    set(handles.checkbox_write,'Value',1)
end

% --- Executes on button press in checkbox_stochastic.
function checkbox_frap_Callback(hObject, eventdata, handles)
global frapch
global ncols
global loadch

if and(ncols==5,and(loadch==1,get(handles.checkbox_frap,'Value')))
    set(handles.checkbox_frap,'Value',0)
    warndlg('FRAP is only supported with a six colum input file.','FRAP Invalid')
    return
end
if frapch==0
    set(handles.checkbox_frap,'Value',0)
    warndlg('First select a FRAP region.','Select Region')
    return
end

% --- Executes on button press in checkbox_log.
function checkbox_log_Callback(hObject, eventdata, handles)
set(handles.checkbox_linear,'Value',0)
set(handles.checkbox_log,'Value',1)
guidata(hObject,handles)

% --- Executes on button press in checkbox_linear.
function checkbox_linear_Callback(hObject, eventdata, handles)
set(handles.checkbox_linear,'Value',1)
set(handles.checkbox_log,'Value',0)
guidata(hObject,handles)

% --- Executes on button press in checkbox_disk.
function checkbox_disk_Callback(hObject, eventdata, handles)
set(handles.checkbox_disk,'Value',1)
set(handles.checkbox_ram,'Value',0)
guidata(hObject,handles)

% --- Executes on button press in checkbox_ram.
function checkbox_ram_Callback(hObject, eventdata, handles)
set(handles.checkbox_disk,'Value',0)
set(handles.checkbox_ram,'Value',1)
guidata(hObject,handles)

% --- Executes on button press in checkbox_histogram.
function checkbox_histogram_Callback(hObject, eventdata, handles)
if get(handles.checkbox_histogram,'Value')
    set(handles.axes2,'Visible','on')
    set(handles.axes2,'XTick',[])
    set(handles.axes2,'YTick',[])
    set(handles.axes2,'Color',[0.9 0.9 0.9])
    set(handles.axes2,'Box','on')
else
    set(handles.axes2,'Visible','off')
end
guidata(hObject,handles)

% --- Executes on button press in checkbox_poiss.
function checkbox_poiss_Callback(hObject, eventdata, handles)
guidata(hObject,handles)

% --- Executes on button press in checkbox_noise.
function checkbox_noise_Callback(hObject, eventdata, handles)
set(handles.checkbox_custom,'Value',0)
guidata(hObject,handles)

% --- Executes on button press in checkbox_custom.
function checkbox_custom_Callback(hObject, eventdata, handles)
global noisech
if noisech==1
    if get(handles.checkbox_custom,'Value')
        set(handles.checkbox_custom,'Value',1)
        set(handles.checkbox_noise,'Value',0)
        guidata(hObject,handles)
    end
else
    warndlg('First load a noise spectrum TIF or text file.','Load Noise File')
    set(handles.checkbox_custom,'Value',0)
    return
end

% --- Executes on button press in checkbox_simZ.
function checkbox_simZ_Callback(hObject, eventdata, handles)
if get(handles.checkbox_simZ,'Value')==0
    set(handles.text49,'Enable','off')
    set(handles.edit_znum,'Enable','off')
    set(handles.text21,'Enable','on')
    set(handles.edit_zslice,'Enable','on')
else
    set(handles.text49,'Enable','on')
    set(handles.edit_znum,'Enable','on')
    set(handles.text21,'Enable','off')
    set(handles.edit_zslice,'Enable','off')
end
guidata(hObject,handles)

% --- Executes on button press in checkbox_fix.
function checkbox_fix_Callback(hObject, eventdata, handles)
global simch
global Ikeep
intlow=str2num(get(handles.edit_lowfix,'String'));
inthigh=str2num(get(handles.edit_highfix,'String'));

if get(handles.checkbox_ram,'Value')
    Imins=min(min(min(Ikeep)));
    Imaxs=max(max(max(Ikeep)));
else
    try
        Imins=2^16-1;
        Imaxs=0;
        for i=1:length(imfinfo('BlurLab_temp.tif'));
            Imins=min([Imins,min(min(imread('BlurLab_temp.tif',i)))]);
            Imaxs=max([Imaxs,max(max(imread('BlurLab_temp.tif',i)))]);
        end
    catch er1
        Imins=0;
        Imaxs=2^16-1;
        warndlg('The temporary stack file is missing, bit values used for contrast.','Missing Temp File')
    end
end

if and(simch==1,get(handles.checkbox_fix,'Value'))
    intlow=round(100*Imins)/100;
    inthigh=round(100*Imaxs)/100;
    
    if get(handles.checkbox_usenoise,'Value')
        %apply scale factor
        sfc=str2num(get(handles.edit_sfc,'String'));
        intlow=intlow*sfc;
        inthigh=inthigh*sfc;
        
        %add Gaussian noise
        if get(handles.checkbox_noise,'Value')
            mu_gauss=str2num(get(handles.edit_mean,'String'));
            std_gauss=str2num(get(handles.edit_std,'String'));
            intlow=intlow+mu_gauss-2.5*std_gauss;
            inthigh=inthigh+mu_gauss+2.5*std_gauss;
        end
        
        %add Poisson/Shot noise
        if get(handles.checkbox_poiss,'Value')
            lam=str2num(get(handles.edit_poiss,'String'));
            basal=str2num(get(handles.edit_basal,'String'));
            intlow=intlow+lam-lam+basal;
            inthigh=inthigh+lam+2.5*lam+basal;
        end
    end
    
    set(handles.edit_lowfix,'String',num2str(intlow));
    set(handles.edit_highfix,'String',num2str(inthigh));
end
edit_slidenum_Callback(hObject, eventdata, handles)
guidata(hObject,handles)

% --- Executes on button press in checkbox_8bit.
function checkbox_8bit_Callback(hObject, eventdata, handles)
set(handles.checkbox_8bit,'Value',1)
set(handles.checkbox_16bit,'Value',0)
guidata(hObject,handles)

% --- Executes on button press in checkbox_16bit.
function checkbox_16bit_Callback(hObject, eventdata, handles)
set(handles.checkbox_8bit,'Value',0)
set(handles.checkbox_16bit,'Value',1)
guidata(hObject,handles)

% --- Executes on button press in checkbox_hot.
function checkbox_hot_Callback(hObject, eventdata, handles)
set(handles.checkbox_hot,'Value',1)
set(handles.checkbox_jet,'Value',0)
set(handles.checkbox_gray,'Value',0)
colormap(hot)
guidata(hObject,handles)

% --- Executes on button press in checkbox_jet.
function checkbox_jet_Callback(hObject, eventdata, handles)
set(handles.checkbox_hot,'Value',0)
set(handles.checkbox_jet,'Value',1)
set(handles.checkbox_gray,'Value',0)
colormap(jet)
guidata(hObject,handles)

% --- Executes on button press in checkbox_gray.
function checkbox_gray_Callback(hObject, eventdata, handles)
set(handles.checkbox_hot,'Value',0)
set(handles.checkbox_jet,'Value',0)
set(handles.checkbox_gray,'Value',1)
colormap(gray)
guidata(hObject,handles)

%**************************************************************************
%***************POP UP MENU************************************************

% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
set(handles.uipanel_output,'Visible','off')
set(handles.uipanel_size,'Visible','off')
set(handles.uipanel_frames,'Visible','off')
set(handles.uipanel_zstacking,'Visible','off')
set(handles.uipanel_randinput,'Visible','off')
set(handles.uipanel_frap,'Visible','off')
set(handles.uipanel_noise,'Visible','off')
set(handles.uipanel_subnoise,'Visible','off')
set(handles.uipanel_tigercreate,'Visible','off')
val1=get(handles.popupmenu1,'Value');
switch val1
    case 1
        set(handles.uipanel_frames,'Visible','on')
    case 2
        set(handles.uipanel_size,'Visible','on')
    case 3
        set(handles.uipanel_noise,'Visible','on')
        set(handles.uipanel_subnoise,'Visible','on')
    case 4
        set(handles.uipanel_zstacking,'Visible','on')
    case 5
        set(handles.uipanel_output,'Visible','on')
    case 6
        set(handles.uipanel_frap,'Visible','on')
    case 7
        set(handles.uipanel_randinput,'Visible','on')
    case 8
        set(handles.uipanel_tigercreate,'Visible','on')
end
% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%************************************



function Tedit_points_Callback(hObject, eventdata, handles)
% hObject    handle to Tedit_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tedit_points as text
%        str2double(get(hObject,'String')) returns contents of Tedit_points as a double


% --- Executes during object creation, after setting all properties.
function Tedit_points_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tedit_points (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tedit_meanI_Callback(hObject, eventdata, handles)
% hObject    handle to Tedit_meanI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tedit_meanI as text
%        str2double(get(hObject,'String')) returns contents of Tedit_meanI as a double


% --- Executes during object creation, after setting all properties.
function Tedit_meanI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tedit_meanI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tedit_D_Callback(hObject, eventdata, handles)
% hObject    handle to Tedit_D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tedit_D as text
%        str2double(get(hObject,'String')) returns contents of Tedit_D as a double


% --- Executes during object creation, after setting all properties.
function Tedit_D_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tedit_D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tedit_Lx_Callback(hObject, eventdata, handles)
% hObject    handle to Tedit_Lx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tedit_Lx as text
%        str2double(get(hObject,'String')) returns contents of Tedit_Lx as a double


% --- Executes during object creation, after setting all properties.
function Tedit_Lx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tedit_Lx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tedit_Ly_Callback(hObject, eventdata, handles)
% hObject    handle to Tedit_Ly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tedit_Ly as text
%        str2double(get(hObject,'String')) returns contents of Tedit_Ly as a double


% --- Executes during object creation, after setting all properties.
function Tedit_Ly_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tedit_Ly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tedit_Lz_Callback(hObject, eventdata, handles)
% hObject    handle to Tedit_Lz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tedit_Lz as text
%        str2double(get(hObject,'String')) returns contents of Tedit_Lz as a double


% --- Executes during object creation, after setting all properties.
function Tedit_Lz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tedit_Lz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Tedit_frames_Callback(hObject, eventdata, handles)
% hObject    handle to Tedit_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tedit_frames as text
%        str2double(get(hObject,'String')) returns contents of Tedit_frames as a double


% --- Executes during object creation, after setting all properties.
function Tedit_frames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tedit_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox33.
function checkbox33_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox33


% --- Executes on button press in checkbox34.
function checkbox34_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox34


% --- Executes on button press in Tig_tog_merging.
function Tig_tog_merging_Callback(hObject, eventdata, handles)
% hObject    handle to Tig_tog_merging (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch get(hObject,'Value')
    case 1
        set(handles.text_MergeD,'Enable','on')
        set(handles.Tig_MergeD,'Enable','on')
        set(handles.text_CMerge,'Enable','on')
        set(handles.Tig_CMerge,'Enable','on')
    case 0
        set(handles.text_MergeD,'Enable','off')
        set(handles.Tig_MergeD,'Enable','off')
        set(handles.text_CMerge,'Enable','off')
        set(handles.Tig_CMerge,'Enable','off')
end


% --- Executes on button press in Tig_tog_spotapp.
function Tig_tog_spotapp_Callback(hObject, eventdata, handles)
% hObject    handle to Tig_tog_spotapp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get(hObject,'Value')
    case 1
        set(handles.text_Apois,'Enable','on')
        set(handles.Tig_Apois,'Enable','on')
    case 0
        set(handles.text_Apois,'Enable','off')
        set(handles.Tig_Apois,'Enable','off')
end


% --- Executes on button press in Tig_tog_blinking.
function Tig_tog_blinking_Callback(hObject, eventdata, handles)
% hObject    handle to Tig_tog_blinking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get(hObject,'Value')
    case 1
        set(handles.text_Nblink,'Enable','on')
        set(handles.Tig_Nblink,'Enable','on')
        set(handles.text_Tblink,'Enable','on')
        set(handles.Tig_Tblink,'Enable','on')
        set(handles.text_Cblink1,'Enable','on')
        set(handles.Tig_Cblink1,'Enable','on')
        set(handles.text_Cblink2,'Enable','on')
        set(handles.Tig_Cblink2,'Enable','on')
    case 0
        set(handles.text_Nblink,'Enable','off')
        set(handles.Tig_Nblink,'Enable','off')
        set(handles.text_Tblink,'Enable','off')
        set(handles.Tig_Tblink,'Enable','off')
        set(handles.text_Cblink1,'Enable','off')
        set(handles.Tig_Cblink1,'Enable','off')
        set(handles.text_Cblink2,'Enable','off')
        set(handles.Tig_Cblink2,'Enable','off')
end

% --- Executes on button press in Tig_tog_splitting.
function Tig_tog_splitting_Callback(hObject, eventdata, handles)
% hObject    handle to Tig_tog_splitting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
switch get(hObject,'Value')
    case 1
        set(handles.text_CSplit,'Enable','on')
        set(handles.Tig_CSplit,'Enable','on')
        set(handles.text_DSplit,'Enable','on')
        set(handles.Tig_DSplit,'Enable','on')
    case 0
        set(handles.text_CSplit,'Enable','off')
        set(handles.Tig_CSplit,'Enable','off')
        set(handles.text_DSplit,'Enable','off')
        set(handles.Tig_DSplit,'Enable','off')
end


% --- Executes on button press in Tig_tog_spotdis.
function Tig_tog_spotdis_Callback(hObject, eventdata, handles)
% hObject    handle to Tig_tog_spotdis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch get(hObject,'Value')
    case 1
        set(handles.text_CDis,'Enable','on')
        set(handles.Tig_CDis,'Enable','on')
    case 0
        set(handles.text_CDis,'Enable','off')
        set(handles.Tig_CDis,'Enable','off')
end



function Tig_SigmaR_Callback(hObject, eventdata, handles)
% hObject    handle to Tig_SigmaR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tig_SigmaR as text
%        str2double(get(hObject,'String')) returns contents of Tig_SigmaR as a double


% --- Executes during object creation, after setting all properties.
function Tig_SigmaR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tig_SigmaR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tig_SigmaA_Callback(hObject, eventdata, handles)
% hObject    handle to Tig_SigmaA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tig_SigmaA as text
%        str2double(get(hObject,'String')) returns contents of Tig_SigmaA as a double


% --- Executes during object creation, after setting all properties.
function Tig_SigmaA_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tig_SigmaA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tig_SigmaI_Callback(hObject, eventdata, handles)
% hObject    handle to Tig_SigmaI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tig_SigmaI as text
%        str2double(get(hObject,'String')) returns contents of Tig_SigmaI as a double


% --- Executes during object creation, after setting all properties.
function Tig_SigmaI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tig_SigmaI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tig_MergeD_Callback(hObject, eventdata, handles)
% hObject    handle to Tig_MergeD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tig_MergeD as text
%        str2double(get(hObject,'String')) returns contents of Tig_MergeD as a double


% --- Executes during object creation, after setting all properties.
function Tig_MergeD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tig_MergeD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tig_CMerge_Callback(hObject, eventdata, handles)
% hObject    handle to Tig_CMerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tig_CMerge as text
%        str2double(get(hObject,'String')) returns contents of Tig_CMerge as a double


% --- Executes during object creation, after setting all properties.
function Tig_CMerge_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tig_CMerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tig_CSplit_Callback(hObject, eventdata, handles)
% hObject    handle to Tig_CSplit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tig_CSplit as text
%        str2double(get(hObject,'String')) returns contents of Tig_CSplit as a double


% --- Executes during object creation, after setting all properties.
function Tig_CSplit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tig_CSplit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tig_DSplit_Callback(hObject, eventdata, handles)
% hObject    handle to Tig_DSplit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tig_DSplit as text
%        str2double(get(hObject,'String')) returns contents of Tig_DSplit as a double


% --- Executes during object creation, after setting all properties.
function Tig_DSplit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tig_DSplit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tig_CDis_Callback(hObject, eventdata, handles)
% hObject    handle to Tig_CDis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tig_CDis as text
%        str2double(get(hObject,'String')) returns contents of Tig_CDis as a double


% --- Executes during object creation, after setting all properties.
function Tig_CDis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tig_CDis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tig_Apois_Callback(hObject, eventdata, handles)
% hObject    handle to Tig_Apois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tig_Apois as text
%        str2double(get(hObject,'String')) returns contents of Tig_Apois as a double


% --- Executes during object creation, after setting all properties.
function Tig_Apois_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tig_Apois (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tig_Nblink_Callback(hObject, eventdata, handles)
% hObject    handle to Tig_Nblink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tig_Nblink as text
%        str2double(get(hObject,'String')) returns contents of Tig_Nblink as a double


% --- Executes during object creation, after setting all properties.
function Tig_Nblink_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tig_Nblink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tig_Tblink_Callback(hObject, eventdata, handles)
% hObject    handle to Tig_Tblink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tig_Tblink as text
%        str2double(get(hObject,'String')) returns contents of Tig_Tblink as a double


% --- Executes during object creation, after setting all properties.
function Tig_Tblink_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tig_Tblink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tig_Cblink1_Callback(hObject, eventdata, handles)
% hObject    handle to Tig_Cblink1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tig_Cblink1 as text
%        str2double(get(hObject,'String')) returns contents of Tig_Cblink1 as a double


% --- Executes during object creation, after setting all properties.
function Tig_Cblink1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tig_Cblink1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Tig_Cblink2_Callback(hObject, eventdata, handles)
% hObject    handle to Tig_Cblink2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tig_Cblink2 as text
%        str2double(get(hObject,'String')) returns contents of Tig_Cblink2 as a double


% --- Executes during object creation, after setting all properties.
function Tig_Cblink2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tig_Cblink2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in Tigercreate_Generate.
function Tigercreate_Generate_Callback(hObject, eventdata, handles)
% hObject    handle to Tigercreate_Generate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ini.SigmaR = str2double(get(handles.Tig_SigmaR,'String'));   % Variance for normal dis. used for change in R (displacement)
ini.SigmaA = str2double(get(handles.Tig_SigmaA,'String'));   % Variance for normal dis. used for change in A (angle)
ini.SigmaI = str2double(get(handles.Tig_SigmaI,'String'));   % Variance for normal dis. used for change in I (Intensity)
ini.MergeD = str2double(get(handles.Tig_MergeD,'String'));   % Distance between two spots needed for merging
ini.CMerge = str2double(get(handles.Tig_CMerge,'String'));   % Chance for a merging event (if within distance)
ini.CSplit = str2double(get(handles.Tig_CSplit,'String'));   % Chance for a splitting event
ini.DSplit = str2double(get(handles.Tig_DSplit,'String'));   % Change in displacement for one particle
ini.CDis = str2double(get(handles.Tig_CDis,'String'));       % Chhance for a disappearence event
ini.Apois = str2double(get(handles.Tig_Apois,'String'));     % Number of spots appearences (poisson distribution)
ini.Nblink = str2double(get(handles.Tig_Nblink,'String'));   % Number of blinking spots
ini.Tblink = str2double(get(handles.Tig_Tblink,'String'));   % Average blinking duration (in frames)
ini.Cblink1 = str2double(get(handles.Tig_Cblink1,'String')); % Chance of becoming a blinking spot (new spots)
ini.Cblink2 = str2double(get(handles.Tig_Cblink2,'String')); % Chance of becoming a stable spot (blinking spots)
ini.MinI = str2double(get(handles.Tig_MinI,'String'));        % Minimal value of intensity

ini.tog_merging = get(handles.Tig_tog_merging,'Value');
ini.tog_split = get(handles.Tig_tog_splitting,'Value');
ini.tog_spotdis = get(handles.Tig_tog_spotdis,'Value');
ini.tog_spotapp = get(handles.Tig_tog_spotapp,'Value');
ini.tog_blinking = get(handles.Tig_tog_blinking,'Value');
ini.tog_dimer = get(handles.Tig_tog_dimer,'Value');
ini.Tigfolder = uigetdir(pwd);

nframes = str2double(get(handles.Tedit_frames,'String'));
pts = str2double(get(handles.Tedit_points,'String'));
meanI = str2double(get(handles.Tedit_meanI,'String'));
D = str2double(get(handles.Tedit_D,'String'));
Lx = str2double(get(handles.Tedit_Lx,'String'));
Ly = str2double(get(handles.Tedit_Ly,'String'));
Lz = str2double(get(handles.Tedit_Lz,'String'));


Tigercreate(nframes,pts,meanI,D,Lx,Ly,Lz,ini)

set(handles.popupmenu1,'Value',1);



function Tig_MinI_Callback(hObject, eventdata, handles)
% hObject    handle to Tig_MinI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tig_MinI as text
%        str2double(get(hObject,'String')) returns contents of Tig_MinI as a double


% --- Executes during object creation, after setting all properties.
function Tig_MinI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tig_MinI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Tig_tog_dimer.
function Tig_tog_dimer_Callback(hObject, eventdata, handles)
% hObject    handle to Tig_tog_dimer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Tig_tog_dimer
