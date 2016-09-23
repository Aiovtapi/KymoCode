function varargout = LionfitcontrolUI(varargin)
% LIONFITCONTROLUI MATLAB code for LionfitcontrolUI.fig
%      LIONFITCONTROLUI, by itself, creates a new LIONFITCONTROLUI or raises the existing
%      singleton*.
%
%      H = LIONFITCONTROLUI returns the handle to a new LIONFITCONTROLUI or the handle to
%      the existing singleton*.
%
%      LIONFITCONTROLUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LIONFITCONTROLUI.M with the given input arguments.
%
%      LIONFITCONTROLUI('Property','Value',...) creates a new LIONFITCONTROLUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LionfitcontrolUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LionfitcontrolUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LionfitcontrolUI

% Last Modified by GUIDE v2.5 31-Aug-2016 18:25:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LionfitcontrolUI_OpeningFcn, ...
                   'gui_OutputFcn',  @LionfitcontrolUI_OutputFcn, ...
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


% --- Executes just before LionfitcontrolUI is made visible.
function LionfitcontrolUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LionfitcontrolUI (see VARARGIN)

handles=guidata(hObject); 

handles.IPTP = varargin{1};
handles.Bacpics = varargin{2};
handles.flimgname = varargin{3};
handles.cells = varargin{4};
handles.chan = varargin{5};
handles.exp = varargin{6};


handles.lioncropindex = 0;

switch handles.exp
    case 'Agar'
        Xout = GaussFitSimedit_AgarControl(handles.Bacpics,1,handles.IPTP,handles.lioncropindex);
    case 'Kymo'
        Xout = GaussFitSimedit_KymoControl(handles.Bacpics,1,handles.IPTP,handles.lioncropindex);
end

title = strcat('Tweak Intensity Peak Threshold for ',handles.flimgname{handles.chan});
set(handles.text_title,'String',title);

set(handles.edit_IPT,'String',num2str(handles.IPTP))
handles.ScaleValue = 0.1;
handles.Cellnumber = 1;

imagesc(handles.Bacpics{1,1},'parent',handles.axes1)
axes(handles.axes1);
hold on
if ~isempty(Xout{1})
    plot(Xout{1}(:,1),Xout{1}(:,2),'kx','LineWidth',2)
end
hold off

imagesc(handles.Bacpics{2,1},'parent',handles.axes2)
axes(handles.axes2);
hold on
if ~isempty(Xout{2})
    plot(Xout{2}(:,1),Xout{2}(:,2),'kx','LineWidth',2)
end
hold off

imagesc(handles.Bacpics{3,1},'parent',handles.axes3)
axes(handles.axes3);
hold on
if ~isempty(Xout{3})
plot(Xout{3}(:,1),Xout{3}(:,2),'kx','LineWidth',2)
end
hold off

% Choose default command line output for LionfitcontrolUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AgarUI wait for user response (see UIRESUME)
uiwait(handles.figure1);




% --- Outputs from this function are returned to the command line.
function varargout = LionfitcontrolUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = str2double(handles.edit_IPT.String);

close(handles.figure1)





function edit_IPT_Callback(hObject, eventdata, handles)
% hObject    handle to edit_IPT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_IPT as text
%        str2double(get(hObject,'String')) returns contents of edit_IPT as a double


% --- Executes during object creation, after setting all properties.
function edit_IPT_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_IPT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Scale_up.
function Scale_up_Callback(hObject, eventdata, handles)

handles.IPTP = handles.IPTP + handles.ScaleValue;
set(handles.edit_IPT,'String',num2str(handles.IPTP));

C1 = handles.Cellnumber;
C2 = C1 + 1; 
C3 = C1 + 2;

switch handles.exp
    case 'Agar'
        Xout = GaussFitSimedit_AgarControl(handles.Bacpics,C1,handles.IPTP,handles.lioncropindex);
    case 'Kymo'
        Xout = GaussFitSimedit_KymoControl(handles.Bacpics,C1,handles.IPTP,handles.lioncropindex);
end

imagesc(handles.Bacpics{C1,1},'parent',handles.axes1)
axes(handles.axes1);
hold on
if ~isempty(Xout{1})
    plot(Xout{1}(:,1),Xout{1}(:,2),'kx','LineWidth',2)
end
hold off

imagesc(handles.Bacpics{C2,1},'parent',handles.axes2)
axes(handles.axes2);
hold on
if ~isempty(Xout{2})
    plot(Xout{2}(:,1),Xout{2}(:,2),'kx','LineWidth',2)
end
hold off

imagesc(handles.Bacpics{C3,1},'parent',handles.axes3)
axes(handles.axes3);
hold on
if ~isempty(Xout{3})
    plot(Xout{3}(:,1),Xout{3}(:,2),'kx','LineWidth',2)
end
hold off


guidata(hObject,handles);


% --- Executes on button press in Scale_down.
function Scale_down_Callback(hObject, eventdata, handles)

handles.IPTP = handles.IPTP - handles.ScaleValue;
set(handles.edit_IPT,'String',num2str(handles.IPTP));

C1 = handles.Cellnumber;
C2 = C1 + 1; 
C3 = C1 + 2;

switch handles.exp
    case 'Agar'
        Xout = GaussFitSimedit_AgarControl(handles.Bacpics,C1,handles.IPTP,handles.lioncropindex);
    case 'Kymo'
        Xout = GaussFitSimedit_KymoControl(handles.Bacpics,C1,handles.IPTP,handles.lioncropindex);
end

imagesc(handles.Bacpics{C1,1},'parent',handles.axes1)
axes(handles.axes1);
hold on
if ~isempty(Xout{1})
    plot(Xout{1}(:,1),Xout{1}(:,2),'kx','LineWidth',2)
end
hold off

imagesc(handles.Bacpics{C2,1},'parent',handles.axes2)
axes(handles.axes2);
hold on
if ~isempty(Xout{2})
    plot(Xout{2}(:,1),Xout{2}(:,2),'kx','LineWidth',2)
end
hold off

imagesc(handles.Bacpics{C3,1},'parent',handles.axes3)
axes(handles.axes3);
hold on
if ~isempty(Xout{3})
    plot(Xout{3}(:,1),Xout{3}(:,2),'kx','LineWidth',2)
end
hold off

guidata(hObject,handles);


% --- Executes on button press in Scale_0001.
function Scale_0001_Callback(hObject, eventdata, handles)
handles.ScaleValue = 0.001;
guidata(hObject,handles);


% --- Executes on button press in Scale_001.
function Scale_001_Callback(hObject, eventdata, handles)
handles.ScaleValue = 0.01;
guidata(hObject,handles);


% --- Executes on button press in Scale_01.
function Scale_01_Callback(hObject, eventdata, handles)
handles.ScaleValue = 0.1;
guidata(hObject,handles);


% --- Executes on button press in Scale_1.
function Scale_1_Callback(hObject, eventdata, handles)
handles.ScaleValue = 1;
guidata(hObject,handles);


% --- Executes on button press in Scale_10.
function Scale_10_Callback(hObject, eventdata, handles)
handles.ScaleValue = 10;
guidata(hObject,handles);


% --- Executes on button press in Previous_set.
function Previous_set_Callback(hObject, eventdata, handles)

handles.Cellnumber = handles.Cellnumber - 3;
C1 = handles.Cellnumber;
C2 = C1 + 1; 
C3 = C1 + 2;

switch handles.exp
    case 'Agar'
        Xout = GaussFitSimedit_AgarControl(handles.Bacpics,C1,handles.IPTP,handles.lioncropindex);
    case 'Kymo'
        Xout = GaussFitSimedit_KymoControl(handles.Bacpics,C1,handles.IPTP,handles.lioncropindex);
end

imagesc(handles.Bacpics{C1,1},'parent',handles.axes1)
axes(handles.axes1);
hold on
if ~isempty(Xout{1})
    plot(Xout{1}(:,1),Xout{1}(:,2),'kx','LineWidth',2)
end
hold off

imagesc(handles.Bacpics{C2,1},'parent',handles.axes2)
axes(handles.axes2);
hold on
if ~isempty(Xout{2})
    plot(Xout{2}(:,1),Xout{2}(:,2),'kx','LineWidth',2)
end
hold off

imagesc(handles.Bacpics{C3,1},'parent',handles.axes3)
axes(handles.axes3);
hold on
if ~isempty(Xout{3})
    plot(Xout{3}(:,1),Xout{3}(:,2),'kx','LineWidth',2)
end
hold off

if C1 == 1;
    set(handles.Previous_set,'Visible','off')
end

if C3 <= handles.cells;
    set(handles.Next_set,'Visible','on')
end

guidata(hObject,handles);


% --- Executes on button press in Next_set.
function Next_set_Callback(hObject, eventdata, handles)

handles.Cellnumber = handles.Cellnumber + 3;
C1 = handles.Cellnumber;
C2 = C1 + 1; 
C3 = C1 + 2;

switch handles.exp
    case 'Agar'
        Xout = GaussFitSimedit_AgarControl(handles.Bacpics,C1,handles.IPTP,handles.lioncropindex);
    case 'Kymo'
        Xout = GaussFitSimedit_KymoControl(handles.Bacpics,C1,handles.IPTP,handles.lioncropindex);
end

imagesc(handles.Bacpics{C1,1},'parent',handles.axes1)
axes(handles.axes1);
hold on
if ~isempty(Xout{1})
    plot(Xout{1}(:,1),Xout{1}(:,2),'kx','LineWidth',2)
end
hold off

imagesc(handles.Bacpics{C2,1},'parent',handles.axes2)
axes(handles.axes2);
hold on
if ~isempty(Xout{2})
    plot(Xout{2}(:,1),Xout{2}(:,2),'kx','LineWidth',2)
end
hold off

imagesc(handles.Bacpics{C3,1},'parent',handles.axes3)
axes(handles.axes3);
hold on
if ~isempty(Xout{3})
    plot(Xout{3}(:,1),Xout{3}(:,2),'kx','LineWidth',2)
end
hold off

if C1 >= 4;
    set(handles.Previous_set,'Visible','on')
end

if C3+3 >= handles.cells;
    set(handles.Next_set,'Visible','off')
end

guidata(hObject,handles);


function SaveIPT_Callback(hObject, eventdata, handles)
uiresume(handles.figure1)
% hObject    handle to SaveIPT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
