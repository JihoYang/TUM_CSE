%===========================================================================%
%|                                                                         |%
%|     Graphical User Interface for pre-processing experimental data       |%
%|                                                                         |%
%|     Model: PLOTTING MODULE                                              |%
%|                                                                         |%
%|     Developer:   Jiho Yang                                              |%
%|                  M.Sc. candidate, Computational Scinece & Engineering   |%
%|                  Technische Universität München, Germany                |%
%|                                                                         |%
%|     Work conducted as a student job (HiWi)                              |%
%|     under Seungyong Oh (Dipl.-Ing)                                      |%
%|     at Lehrstuhl für Fördertechnik Materialfluss Logistik               |%
%|     Dept of Mechanical Engineering                                      |%
%|     Technische Universität München, Germany                             |%
%|                                                                         |%
%|     Final update: 9th April 2016                                        |%  
%|                                                                         |%
%===========================================================================%


%Builds basic MATLAB GUI environment. DO NOT MAKE CHANGES
%==================================================================================%
function varargout = faml_plot(varargin)
% faml_plot MATLAB code for faml_plot.fig
%      FAML_PLOT, by itself, creates a new FAML_PLOT or raises the existing
%      singleton*.
%
%      H = FAML_PLOT returns the handle to a new FAML_PLOT or the handle to
%      the existing singleton*.
%
%      FAML_PLOT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FAML_PLOT.M with the given input arguments.
%
%      FAML_PLOT('Property','Value',...) creates a new faml_plot or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before faml_plot_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to faml_plot_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help faml_plot

% Last Modified by GUIDE v2.5 09-Apr-2016 00:02:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @faml_plot_OpeningFcn, ...
                   'gui_OutputFcn',  @faml_plot_OutputFcn, ...
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

% --- Executes just before faml_plot is made visible.
function faml_plot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to faml_plot (see VARARGIN)

% Choose default command line output for faml_plot
handles.output = hObject;

%Displays FML logo
axes(handles.logo);
imshow('Logo.png');

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes faml_plot wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = faml_plot_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%FOLLOWING FUNCTIONS MAY BE MODIFIED
%==================================================================================%

% --- Executes on button press in browse.
function browse_Callback(hObject, eventdata, handles)
% hObject    handle to browse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Bring the selected text file (.xlsx) and save its file and path names
[filename_exp pathname_exp] = uigetfile({'.xlsx'},'File Selector');
fullpathname_exp = strcat(pathname_exp, filename_exp);

%Build full file name and read/import it to MATLAB GUI workspace
loaddata_exp = fullfile(pathname_exp, filename_exp);
data_numeric = importdata(loaddata_exp);

%Save the imported data
handles.data_numeric = data_numeric;
guidata(hObject, handles);

%Plot data in GUI environment
plot(handles.plot, data_numeric(:,1), data_numeric(:,2));
box(handles.plot, 'on');
set(handles.plot, 'FontSize', 15);
set(handles.plot, 'FontName', 'Arial');

%Display full file path directory in GUI environment
set(handles.filepath, 'String', fullpathname_exp);

% --- Executes on button press in configure_plot.
function configure_plot_Callback(hObject, eventdata, handles)
% hObject    handle to configure_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Bring the user defined figure configurations
t = get(handles.title, 'String');
x = get(handles.xlabel, 'String');
y = get(handles.ylabel, 'String');
font = get(handles.font, 'String');
fontsize = str2double(get(handles.fontsize, 'String'));

%Configure the plot in GUI environment
title(handles.plot, t);
xlabel(handles.plot, x);
ylabel(handles.plot, y);
set(handles.plot, 'FontSize', fontsize);
set(handles.plot, 'FontName', font);

%Configure and plot in additional figure
data_numeric = handles.data_numeric;
figure (1)
plot(data_numeric(:,1), data_numeric(:,2));
title(t);
xlabel(x);
ylabel(y);
box on;
set(gca, 'FontSize', fontsize);
set(gca, 'FontName', strcat(font));



%DEFINE GUI COMPONENTS: DO NOT MAKE CHANGES
%==================================================================================%
function title_Callback(hObject, eventdata, handles)
% hObject    handle to title (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of title as text
%        str2double(get(hObject,'String')) returns contents of title as a double

% --- Executes during object creation, after setting all properties.
function title_CreateFcn(hObject, eventdata, handles)
% hObject    handle to title (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function font_Callback(hObject, eventdata, handles)
% hObject    handle to font (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of font as text
%        str2double(get(hObject,'String')) returns contents of font as a double

% --- Executes during object creation, after setting all properties.
function font_CreateFcn(hObject, eventdata, handles)
% hObject    handle to font (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function xlabel_Callback(hObject, eventdata, handles)
% hObject    handle to xlabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of xlabel as text
%        str2double(get(hObject,'String')) returns contents of xlabel as a double

% --- Executes during object creation, after setting all properties.
function xlabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xlabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ylabel_Callback(hObject, eventdata, handles)
% hObject    handle to ylabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of ylabel as text
%        str2double(get(hObject,'String')) returns contents of ylabel as a double

% --- Executes during object creation, after setting all properties.
function ylabel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ylabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function fontsize_Callback(hObject, eventdata, handles)
% hObject    handle to fontsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of fontsize as text
%        str2double(get(hObject,'String')) returns contents of fontsize as a double

% --- Executes during object creation, after setting all properties.
function fontsize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fontsize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%==================================================================================%
