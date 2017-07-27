
%===========================================================================%
%|                                                                         |%
%|     Graphical User Interface for post-processing experimental data      |%
%|                                                                         |%
%|     Model: HYDROPULSER                                                  |%
%|                                                                         |%
%|     Developer:   Jiho Yang (MEng)                                       |%
%|                  M.Sc. candidate, Computational Scinece & Engineering   |%
%|                  Technische Universitat Munchen, Germany                |%
%|                                                                         |%
%|     Work conducted as a student job (HiWi)                              |%
%|     under Seungyong Oh (Dipl.-Ing)                                      |%
%|     at Lehrstuhl fur Fordertechnik Materialfluss Logistik               |%
%|     Dept of Mechanical Engineering                                      |%
%|     Technische Universitat Munchen, Germany                             |%
%|                                                                         |%
%|     Final update: 27th May 2016                                         |%  
%|                                                                         |%
%===========================================================================%


%Builds basic MATLAB GUI environment. DO NOT MAKE CHANGES
%==================================================================================%
function varargout = faml_Hydropulser_GUI(varargin)
% faml_Hydropulser_GUI MATLAB code for faml_faml_Hydropulser_GUI.fig
%      faml_Hydropulser_GUI, by itself, creates a new faml_Hydropulser_GUI
%      or raises the existing singleton*.
%
%      H = faml_Hydropulser_GUI returns the handle to a new
%      faml_Hydropulser_GUI or the handle to the existing singleton*.
%
%      faml_Hydropulser_GUI('CALLBACK',hObject,eventData,handles,...) calls
%      the local function named CALLBACK in faml_Hydropulser_GUI.M with the
%      given input arguments.
%
%      faml_Hydropulser_GUI('Property','Value',...) creates a new
%      faml_Hydropulser_GUI or raises the existing singleton*.  Starting
%      from the left, property value pairs are applied to the GUI before
%      faml_Hydropulser_GUI_OpeningFcn gets called.  An unrecognized
%      property name or invalid value makes property application stop.  All
%      inputs are passed to faml_Hydropulser_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help faml_Hydropulser_GUI

% Last Modified by GUIDE v2.5 27-May-2016 01:06:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @faml_Hydropulser_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @faml_Hydropulser_GUI_OutputFcn, ...
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

% --- Executes just before faml_Hydropulser_GUI is made visible.
function faml_Hydropulser_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn. hObject    handle to
% figure eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA) varargin
% command line arguments to faml_Hydropulser_GUI (see VARARGIN)

% Choose default command line output for faml_Hydropulser_GUI
handles.output = hObject;

%Displays FML logo
axes(handles.logo);
imshow('Logo.png');

%Define simulation data index for multiple simulation data comparison
sim_index = 1;
handles.sim_index = sim_index;
guidata(hObject, handles);

% Update handles structure
% UIWAIT makes faml_Hydropulser_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = faml_Hydropulser_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT); hObject
% handle to figure eventdata  reserved - to be defined in a future version
% of MATLAB handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
%==================================================================================%




%FOLLOWING FUNCTIONS MAY BE MODIFIED
%==================================================================================%

% --- Executes on button press in Browse_ExperimentalData.
% --- Loads and processes experimental data in text file format (.txt)
function Browse_ExperimentalData_Callback(hObject, eventdata, handles)
% hObject    handle to Browse_ExperimentalData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

     %Bring the selected text file (.txt.) and save its file and path names
     [filename pathname] = uigetfile({'.txt'},'File Selector');
     fullpathname = strcat(pathname, filename);
     
     %Build full file name and read/import it to MATLAB GUI workspace
     loaddata = fullfile(pathname, filename);
     data_numeric = dlmread(loaddata,'\t');
     
     %============================================================================%
     %The following lines are used for importing text files with ',' as
     %decimal place indicator
     
     %Reads ',' as strings and converts them into '.'
     
     %data_numeric = [];
 
     %   for i=1:size(data_string,1)
         
     %       data_numeric(i,:) = str2num(char(strrep(data_string(i,:),',','.')'));
          
     %   end
     %============================================================================%
     
     %Define parameters for data plotting
     t_min = min(data_numeric(:,1));
     t_max = max(data_numeric(:,1));
     
     output = {'Force', 'Displacement', 'Temperature'};
     unit   = {' (kN)', ' (mm)', ' (\circC)'};
     
     %Create handles for sharing data between callbacks
     handles.data_exp = data_numeric;
     guidata(hObject, handles);
        
     
     %Plot data in GUI environment
     cla(handles.axes1); cla(handles.axes2); cla(handles.axes3);
     faml_plot_GUI('', data_numeric, output, unit, t_min, t_max, 'Experiment', handles);
 
     %Plot data in separate figure windows
     cla(figure(1)); cla(figure(2)); cla(figure(3));
     faml_plot_figure('', data_numeric, t_min, t_max, output, unit, 'Experiment');
     
     %Display full file path directory in GUI environment
     set(handles.experiment_filepath, 'String', fullpathname);

     
% --- Executes on button press in Browse_SimulationData.
% --- Loads and processes simulation data in text file format (.tab)
function Browse_SimulationData_Callback(hObject, eventdata, handles)
% hObject    handle to Browse_SimulationData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %Bring sim_index. Initial value defined in FAML_HYDROPULSER_GUI_OPENINGFCN
    sim_index = handles.sim_index;
    
    %Bring the selected text file (.dat) and save its file and path names
    [filename_sim pathname_sim] = uigetfile({'.dat'},'File Selector');
    fullpathname_sim = strcat(pathname_sim, filename_sim);

    %Build full file name and read/import it to MATLAB GUI workspace
    loaddata_sim = fullfile(pathname_sim, filename_sim);
    data_cellarray = importdata(loaddata_sim);
    
    %Dynamic allocation of imported simulation data into cell array
    data_numeric_sim{sim_index} = data_cellarray.data;
    
    %Define parameters for data plotting
    t_min_sim{sim_index} = min(data_numeric_sim{sim_index}(:,1));
    t_max_sim{sim_index} = max(data_numeric_sim{sim_index}(:,1));
 
    output_sim = {'Force', 'Displacement'};
    unit_sim   = {' (kN)', ' (mm)'};
    
    %Define matrices & vectors needed for calibration
    time_force{sim_index} = data_numeric_sim{sim_index}(:,1);
    time_disp{sim_index} = data_numeric_sim{sim_index}(:,1);
    disp_force{sim_index} = [data_numeric_sim{sim_index}(:,2), data_numeric_sim{sim_index}(:,3)];
    
    %Create a cell array for saving calibration values
    net_calibration{1,sim_index} = zeros(1,3);
    net_calibration{1,sim_index}(1,3) = 1;
    
    net_calibration{2,sim_index} = zeros(1,3);
    net_calibration{2,sim_index}(1,3) = 1;
    
    net_calibration{3,sim_index} = zeros(1,3);
    net_calibration{3,sim_index}(1,3) = 1;
    
    %Create a cell array for showing the net calibration value
    legend_m_force = cell(2,1);
    legend_m_force(1) = {'Experiment'};
    legend_m_disp = cell(2,1);
    legend_m_disp(1) = {'Experiment'};
    legend_m_disp_force = cell(2,1);
    legend_m_disp_force(1) = {'Experiment'};
     
    %Create handles for sharing data
    handles.net_calibration{1,sim_index} = net_calibration{1,sim_index};
    guidata(hObject, handles);
    
    handles.net_calibration{2,sim_index} = net_calibration{2,sim_index};
    guidata(hObject, handles);
    
    handles.net_calibration{3,sim_index} = net_calibration{3,sim_index};
    guidata(hObject, handles);
    
    handles.legend_m_force{1,1} = legend_m_force{1,1};
    guidata(hObject, handles);
    
    handles.legend_m_disp{1,1} = legend_m_disp{1,1};
    guidata(hObject, handles);
    
    handles.legend_m_disp_force{1,1} = legend_m_disp_force{1,1};
    guidata(hObject, handles);
    
    handles.data_sim{sim_index} = data_numeric_sim{sim_index};
    guidata(hObject, handles);
    
    handles.time_force{sim_index} = time_force{sim_index};
    guidata(hObject, handles);
    
    handles.time_disp{sim_index} = time_disp{sim_index};
    guidata(hObject, handles);
    
    handles.disp_force{sim_index} = disp_force{sim_index};
    guidata(hObject, handles);
    
    handles.t_min_sim{sim_index}=t_min_sim{sim_index};
    guidata(hObject, handles);
      
    handles.t_max_sim{sim_index}=t_max_sim{sim_index};
    guidata(hObject, handles);
    
    %Create and save a matrix for plot legends     
    handles.legend_m_sim{sim_index,1} = strcat('T', num2str(sim_index),' ', 'Force', ' vs Time');
    handles.legend_m_sim{sim_index,2} = strcat('T', num2str(sim_index),' ', 'Displacement', ' vs Time');
    handles.legend_m_sim{sim_index,3} = strcat('T', num2str(sim_index),' ', 'Displacement vs Force');    
    handles.legend_m_sim{sim_index,4} = strcat('T', num2str(sim_index));
    handles.legend_m_sim{sim_index,5} = strcat('T', num2str(sim_index));
    guidata(hObject, handles);
    
    legend_m_sim = handles.legend_m_sim;
    
    
    %Initialise legend_m_defl
    handles.legend_m_force{sim_index+1} = strcat('T',num2str(sim_index), ...
                                                     '(x =  0', ...
                                                     ' y =  0', ...
                                                     ' scale = 1)');
    guidata(hObject, handles);
                                                 
                             
    handles.legend_m_disp{sim_index+1} = strcat('T',num2str(sim_index), ...
                                                     '(x =  0', ...
                                                     ' y =  0', ...
                                                     ' scale = 1)');
    guidata(hObject, handles);
    
    handles.legend_m_disp_force{sim_index+1} = strcat('T',num2str(sim_index), ...
                                                     '(x =  0', ...
                                                     ' y =  0', ...
                                                     ' scale = 1)');
  
    guidata(hObject, handles);
    
    %Plot data in GUI environment
    %cla(handles.axes4); cla(handles.axes5); cla(handles.axes6);
    %cla(handles.axes7); cla(handles.axes8);
    faml_plot_GUI(legend_m_sim, data_numeric_sim{sim_index}, output_sim, unit_sim, t_min_sim{sim_index}, t_max_sim{sim_index}, 'Simulation', handles);
 
    %Plot data in separate figure windows
    %cla(figure(4)); cla(figure(5)); cla(figure(6));
    %cla(figure(7)); cla(figure(8));
    faml_plot_figure(legend_m_sim, data_numeric_sim{sim_index}, t_min_sim{sim_index}, t_max_sim{sim_index}, output_sim, unit_sim, 'Simulation');
    
    %Update simulation data index
    handles.sim_index = sim_index + 1;
    guidata(hObject, handles); 
    
    %Save the imported simulation data in MATLAB base workspace
    handles.data_sim_original{sim_index} = handles.data_sim{sim_index};
    guidata(hObject, handles);
    
    handles.time_force_original{sim_index} = time_force{sim_index};
    guidata(hObject, handles);
    
    handles.time_disp_original{sim_index} = time_disp{sim_index};
    guidata(hObject, handles);
    
    handles.disp_force_original{sim_index} = disp_force{sim_index};
    guidata(hObject, handles);
    
    %Initialise calibrated simulation data cell array for calibration
    data_force_calibrated{sim_index} = [handles.time_force{sim_index}(:), handles.data_sim{sim_index}(:,2)];
    data_disp_calibrated{sim_index} = [handles.time_disp{sim_index}(:), handles.data_sim{sim_index}(:,3)];
    data_disp_force_calibrated{sim_index} = handles.disp_force{sim_index};
    
    %Save initialised calibrated simulation data
    handles.data_force_calibrated{sim_index} = data_force_calibrated{sim_index};
    guidata(hObject, handles);
    handles.data_disp_calibrated{sim_index} = data_disp_calibrated{sim_index};
    guidata(hObject, handles);
    handles.data_disp_force_calibrated{sim_index} = data_disp_force_calibrated{sim_index};
    guidata(hObject, handles);
       
    %Display full file path directory in GUI environment
    set(handles.simulation_filepath, 'String', fullpathname_sim);
  
    
% --- Calibrates simulation data
% --- Calibration made with positive signs (+) (i.e. calibrate value of -1
%     will either subtract the data by -1 or change the sign)
% --- PLEASE MAKE SURE TO SELECT ONLY ONE OPERATION (may give wrong plots)
% --- Executes on button press in calibrate_button.
function calibrate_button_Callback(hObject, eventdata, handles)
% hObject    handle to calibrate_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    %calibrate_value : User specified calibration value from calibration_value handle
    %data_original   : Imported experimental data from Browse_ExperimenetalData callback
    %data_sim        : Imported simulation data from Browse_SimulationData callback
    
    %Bring user specified calibration value from calibration_value callback
    calibrate_value = str2double(get(handles.calibration_value, 'String'));

    %Bring imported experimental data
    data_exp = handles.data_exp;
    
    %Define strings for data plotting
    ptitle = {'Force vs time (calibrated)', 'Displacement vs time (calibrated)', ...
               'Displacement vs Force (calibrated)'}; 
    xoutput = {'Time (sec)', 'Time (sec)', 'Force (kN)'};
    youtput = {'Force (kN)', 'Displacement (mm)', 'Displacement (mm)'};
    
    %Specify which simulation data to calibrate
    if get(handles.T1, 'Value')
        
        data_id = 1;
        
    end
    
    if get(handles.T2, 'Value')
        
       data_id = 2;
        
    end
    
    if get(handles.T3, 'Value')
        
      data_id = 3;
        
    end
    
    if get(handles.T4, 'Value')
        
      data_id = 4;
        
    end
    
    %Bring computed minimum and maximum time from the simulation dataset
    t_min_sim{data_id} = handles.t_min_sim{data_id};
    t_max_sim{data_id} = handles.t_max_sim{data_id};
    
    
    %Calibrate force vs time
    if get(handles.calibrate_force, 'Value')
       
       %Condition statements for shift/scale operations
    
        %Shift in Y axes
        if get(handles.calibration_checkbox_shift, 'Value')
           
           handles.data_sim{data_id}(:,2) = handles.data_sim{data_id}(:,2) + calibrate_value;
           guidata(hObject, handles);
           handles.data_force_calibrated{data_id} = [handles.time_force{data_id}(:), handles.data_sim{data_id}(:,2)];
           guidata(hObject, handles);
           
           handles.net_calibration{1,data_id}(1,2) = handles.net_calibration{1,data_id}(1,2) + calibrate_value;
           guidata(hObject, handles);
                                                   
           handles.legend_m_force{data_id+1} = strcat('T',num2str(data_id), ...
                                                      '(x =   ', num2str(handles.net_calibration{1,data_id}(1,1)), ...
                                                      ' y =   ', num2str(handles.net_calibration{1,data_id}(1,2)), ...
                                                      ' scale =   ', num2str(handles.net_calibration{1,data_id}(1,3)), ')');
                                                                   
           guidata(hObject, handles);
 
           faml_plot_figure_calibration(handles.legend_m_force, data_exp, handles.data_force_calibrated, ptitle, xoutput, youtput, t_min_sim{data_id}, t_max_sim{data_id}, handles);  
           faml_plot_GUI_calibration(handles.legend_m_force, data_exp, handles.data_force_calibrated, ptitle, xoutput, youtput, t_min_sim{data_id}, t_max_sim{data_id}, handles);
    
        end
        
        %Shift in X axes
        if get(handles.calibration_checkbox_shift_time, 'Value')
           
           handles.time_force{data_id} = handles.time_force{data_id} + calibrate_value;
           guidata(hObject, handles);
           handles.data_force_calibrated{data_id} = [handles.time_force{data_id}, handles.data_sim{data_id}(:,2)];
           guidata(hObject, handles);
           
           handles.net_calibration{1,data_id}(1,1) = handles.net_calibration{1,data_id}(1,1) + calibrate_value;
           guidata(hObject, handles);
           
           handles.legend_m_force{data_id+1} = strcat('T',num2str(data_id), ...
                                                      '(x =   ', num2str(handles.net_calibration{1,data_id}(1,1)), ...
                                                      ' y =   ', num2str(handles.net_calibration{1,data_id}(1,2)), ...
                                                      ' scale =   ', num2str(handles.net_calibration{1,data_id}(1,3)), ')');
                                                                   
           guidata(hObject, handles);
           
           t_min_shifted{data_id} = min(handles.data_force_calibrated{data_id}(:,1));
           t_max_shifted{data_id} = max(handles.data_force_calibrated{data_id}(:,1));
           
           t_domain_min{data_id} = [t_min_sim{data_id}, t_min_shifted{data_id}];
           t_domain_max{data_id} = [t_max_sim{data_id}, t_max_shifted{data_id}];
           
           t_min_calibrated{data_id} = t_domain_min{data_id}(t_domain_min{data_id} >= t_min_sim{data_id});
           t_max_calibrated{data_id} = t_domain_max{data_id}(t_domain_max{data_id} <= t_max_sim{data_id});
           
                %IF statements for constraining x-axes
                if numel(t_min_calibrated{data_id}) == 2 
                    
                   t_min_calibrated{data_id} = t_min_shifted{data_id};
                   
                end
                
                if numel(t_max_calibrated{data_id}) == 2
                    
                   t_max_calibrated{data_id} = t_max_shifted{data_id};
                   
                end
                
                if t_min_calibrated{data_id} > t_max_sim{data_id} || t_max_calibrated{data_id} < t_min_sim{data_id}
                    
                   display('Calibrated data is out of domain, please reset');
                   return;
                    
                end
            
           faml_plot_figure_calibration(handles.legend_m_force, data_exp, handles.data_force_calibrated, ptitle, xoutput, youtput, t_min_calibrated{data_id}, t_max_calibrated{data_id}, handles);  
           faml_plot_GUI_calibration(handles.legend_m_force, data_exp, handles.data_force_calibrated, ptitle, xoutput, youtput, t_min_sim{data_id}, t_max_sim{data_id}, handles);
    
        end
        
        %Scale in Y axes
        if get(handles.calibration_checkbox_scale, 'Value')
            
           handles.data_sim{data_id}(:,2) = handles.data_sim{data_id}(:,2)*calibrate_value;
           guidata(hObject, handles);
           handles.data_force_calibrated{data_id} = [handles.time_force{data_id}, handles.data_sim{data_id}(:,2)];
           guidata(hObject, handles);
           
           handles.net_calibration{1,data_id}(1,3) = handles.net_calibration{1,data_id}(1,3) * calibrate_value;
           guidata(hObject, handles);
           
           handles.legend_m_force{data_id+1} = strcat('T',num2str(data_id), ...
                                                      '(x =   ', num2str(handles.net_calibration{1,data_id}(1,1)), ...
                                                      ' y =   ', num2str(handles.net_calibration{1,data_id}(1,2)), ...
                                                      ' scale =   ', num2str(handles.net_calibration{1,data_id}(1,3)), ')');
                                                                   
           guidata(hObject, handles);
      
           faml_plot_figure_calibration(handles.legend_m_force, data_exp, handles.data_force_calibrated, ptitle, xoutput, youtput, t_min_sim{data_id}, t_max_sim{data_id}, handles);  
           faml_plot_GUI_calibration(handles.legend_m_force, data_exp, handles.data_force_calibrated, ptitle, xoutput, youtput, t_min_sim{data_id}, t_max_sim{data_id}, handles);
            
        end
        
        %Display warning messages when multiple operations are selected
        if get(handles.calibration_checkbox_shift, 'Value') && ...
           get(handles.calibration_checkbox_scale, 'Value') || ...
           get(handles.calibration_checkbox_shift, 'Value') && ...
           get(handles.calibration_checkbox_shift_time, 'Value') || ...
           get(handles.calibration_checkbox_shift_time, 'Value') && ...
           get(handles.calibration_checkbox_scale, 'Value')
            
           display('Multiple operation selected. Please choose either Scale or Shift');
           return;
           
        end
        
        if get(handles.calibrate_force, 'Value') && ...
           get(handles.calibrate_disp, 'Value') || ...
           get(handles.calibrate_force, 'Value') && ...
           get(handles.calibrate_disp_force, 'Value') || ...
           get(handles.calibrate_disp, 'Value') && ...
           get(handles.calibrate_disp_force, 'Value')
            
           display('Multiple variable selected. Please choose one variable only');
           return;
     
        end
        
    end
  
    %Calibrate Displacement vs Time
    if get(handles.calibrate_disp, 'Value')
       
       %Condition statements for shift/scale operations
    
        %Shift in Y axes
        if get(handles.calibration_checkbox_shift, 'Value')
           
           handles.data_sim{data_id}(:,3) = handles.data_sim{data_id}(:,3) + calibrate_value;
           guidata(hObject, handles);
           handles.data_disp_calibrated{data_id} = [handles.time_disp{data_id}, handles.data_sim{data_id}(:,3)];
           guidata(hObject, handles);
           
           handles.net_calibration{2,data_id}(1,2) = handles.net_calibration{1,data_id}(1,2) + calibrate_value;
           guidata(hObject, handles);
           
           handles.legend_m_disp{data_id+1} = strcat('T',num2str(data_id), ...
                                                     '(x =   ', num2str(handles.net_calibration{2,data_id}(1,1)), ...
                                                     ' y =   ', num2str(handles.net_calibration{2,data_id}(1,2)), ...
                                                     ' scale =   ', num2str(handles.net_calibration{2,data_id}(1,3)), ')');
                                                                   
           guidata(hObject, handles);
            
           faml_plot_figure_calibration(handles.legend_m_disp, data_exp, handles.data_disp_calibrated, ptitle, xoutput, youtput, t_min_sim{data_id}, t_max_sim{data_id}, handles);  
           faml_plot_GUI_calibration(handles.legend_m_disp, data_exp, handles.data_disp_calibrated, ptitle, xoutput, youtput, t_min_sim{data_id}, t_max_sim{data_id}, handles);
    
        end
        
        %Shift in X axes
        if get(handles.calibration_checkbox_shift_time, 'Value')
           
           handles.time_disp{data_id} = handles.time_disp{data_id} + calibrate_value;
           guidata(hObject, handles);
           handles.data_disp_calibrated{data_id} = [handles.time_disp{data_id}, handles.data_sim{data_id}(:,3)];
           guidata(hObject, handles);
           
           handles.net_calibration{2,data_id}(1,1) = handles.net_calibration{2,data_id}(1,1) + calibrate_value;
           guidata(hObject, handles);
           
           handles.legend_m_disp{data_id+1} = strcat('T',num2str(data_id), ...
                                                     '(x =   ', num2str(handles.net_calibration{2,data_id}(1,1)), ...
                                                     ' y =   ', num2str(handles.net_calibration{2,data_id}(1,2)), ...
                                                     ' scale =   ', num2str(handles.net_calibration{2,data_id}(1,3)), ')');
                                                                   
           guidata(hObject, handles);
           
           t_min_shifted{data_id} = min(handles.data_disp_calibrated{data_id}(:,1));
           t_max_shifted{data_id} = max(handles.data_disp_calibrated{data_id}(:,1));
           
           t_domain_min{data_id} = [t_min_sim{data_id}, t_min_shifted{data_id}];
           t_domain_max{data_id} = [t_max_sim{data_id}, t_max_shifted{data_id}];
           
           t_min_calibrated{data_id} = t_domain_min{data_id}(t_domain_min{data_id} >= t_min_sim{data_id});
           t_max_calibrated{data_id} = t_domain_max{data_id}(t_domain_max{data_id} <= t_max_sim{data_id});
           
                %IF statements for constraining x-axes
                if numel(t_min_calibrated{data_id}) == 2 
                    
                   t_min_calibrated{data_id} = t_min_shifted{data_id};
                   
                end
                
                if numel(t_max_calibrated{data_id}) == 2
                    
                   t_max_calibrated{data_id} = t_max_shifted{data_id};
                   
                end
                
                if t_min_calibrated{data_id} > t_max_sim{data_id} || t_max_calibrated{data_id} < t_min_sim{data_id}
                    
                   display('Calibrated data is out of domain, please reset');
                   return;
                    
                end
            
           faml_plot_figure_calibration(handles.legend_m_disp, data_exp, handles.data_disp_calibrated, ptitle, xoutput, youtput, t_min_calibrated{data_id}, t_max_calibrated{data_id}, handles);  
           faml_plot_GUI_calibration(handles.legend_m_disp, data_exp, handles.data_disp_calibrated, ptitle, xoutput, youtput, t_min_sim{data_id}, t_max_sim{data_id}, handles);
    
        end
        
        %Scale in Y axes
        if get(handles.calibration_checkbox_scale, 'Value')
            
           handles.data_sim{data_id}(:,3) = handles.data_sim{data_id}(:,3)*calibrate_value;
           guidata(hObject, handles);
           handles.data_disp_calibrated{data_id} = [handles.time_disp{data_id}, handles.data_sim{data_id}(:,3)];
           guidata(hObject, handles);
           
           handles.net_calibration{2,data_id}(1,3) = handles.net_calibration{2,data_id}(1,3) * calibrate_value;
           guidata(hObject, handles);
           
           handles.legend_m_disp{data_id+1} = strcat('T',num2str(data_id), ...
                                                     '(x =   ', num2str(handles.net_calibration{2,data_id}(1,1)), ...
                                                     ' y =   ', num2str(handles.net_calibration{2,data_id}(1,2)), ...
                                                     ' scale =   ', num2str(handles.net_calibration{2,data_id}(1,3)), ')');
                                                                   
           guidata(hObject, handles);
      
           faml_plot_figure_calibration(handles.legend_m_disp, data_exp, handles.data_disp_calibrated, ptitle, xoutput, youtput, t_min_sim{data_id}, t_max_sim{data_id}, handles);  
           faml_plot_GUI_calibration(handles.legend_m_disp, data_exp, handles.data_disp_calibrated, ptitle, xoutput, youtput, t_min_sim{data_id}, t_max_sim{data_id}, handles);
            
        end
        
        %Display warning messages when multiple operations are selected
        if get(handles.calibration_checkbox_shift, 'Value') && ...
           get(handles.calibration_checkbox_scale, 'Value') || ...
           get(handles.calibration_checkbox_shift, 'Value') && ...
           get(handles.calibration_checkbox_shift_time, 'Value') || ...
           get(handles.calibration_checkbox_shift_time, 'Value') && ...
           get(handles.calibration_checkbox_scale, 'Value')
            
           display('Multiple operation selected. Please choose either Scale or Shift');
           return;
           
        end
        
        if get(handles.calibrate_force, 'Value') && ...
           get(handles.calibrate_disp, 'Value') || ...
           get(handles.calibrate_force, 'Value') && ...
           get(handles.calibrate_disp_force, 'Value') || ...
           get(handles.calibrate_disp, 'Value') && ...
           get(handles.calibrate_disp_force, 'Value')
            
           display('Multiple variable selected. Please choose one variable only');
           return;
     
        end
     
    end
    
    %Calibrate Displacement vs Force
    if get(handles.calibrate_disp_force, 'Value')
       
        %Condition statements for shift/scale operations
    
        %Shift in Y axes
        if get(handles.calibration_checkbox_shift, 'Value')
           
           handles.disp_force{data_id}(:,2) = handles.disp_force{data_id}(:,2) + calibrate_value;
           guidata(hObject, handles);
           handles.data_disp_force_calibrated{data_id} = handles.disp_force{data_id};
           guidata(hObject, handles);
           
           handles.net_calibration{3,data_id}(1,2) = handles.net_calibration{3,data_id}(1,2) + calibrate_value;
           guidata(hObject, handles);
           
           handles.legend_m_disp_force{data_id+1} = strcat('T',num2str(data_id), ...
                                                           '(x =   ', num2str(handles.net_calibration{3,data_id}(1,1)), ...
                                                           ' y =   ', num2str(handles.net_calibration{3,data_id}(1,2)), ...
                                                           ' scale =   ', num2str(handles.net_calibration{3,data_id}(1,3)), ')');
                                                                   
           guidata(hObject, handles);
            
           faml_plot_figure_calibration(handles.legend_m_disp_force, data_exp, handles.data_disp_force_calibrated, ptitle, xoutput, youtput, t_min_sim{data_id}, t_max_sim{data_id}, handles);  
           faml_plot_GUI_calibration(handles.legend_m_disp_force, data_exp, handles.data_disp_force_calibrated, ptitle, xoutput, youtput, t_min_sim{data_id}, t_max_sim{data_id}, handles);
    
        end
        
        %Shift in X axes
        if get(handles.calibration_checkbox_shift_time, 'Value')
           
           handles.disp_force{data_id}(:,1) = handles.disp_force{data_id}(:,1) + calibrate_value;
           guidata(hObject, handles);
           handles.data_disp_force_calibrated{data_id} = handles.disp_force{data_id};
           guidata(hObject, handles);
           
           handles.net_calibration{3,data_id}(1,1) = handles.net_calibration{3,data_id}(1,1) + calibrate_value;
           guidata(hObject, handles);
           
           handles.legend_m_disp_force{data_id+1} = strcat('T',num2str(data_id), ...
                                                           '(x =   ', num2str(handles.net_calibration{3,data_id}(1,1)), ...
                                                           ' y =   ', num2str(handles.net_calibration{3,data_id}(1,2)), ...
                                                           ' scale =   ', num2str(handles.net_calibration{3,data_id}(1,3)), ')');
                                                                   
           guidata(hObject, handles);
           
           t_min_shifted{data_id} = min(handles.data_disp_force_calibrated{data_id}(:,1));
           t_max_shifted{data_id} = max(handles.data_disp_force_calibrated{data_id}(:,1));
           
           t_domain_min{data_id} = [t_min_sim{data_id}, t_min_shifted{data_id}];
           t_domain_max{data_id} = [t_max_sim{data_id}, t_max_shifted{data_id}];
           
           t_min_calibrated{data_id} = t_domain_min{data_id}(t_domain_min{data_id} >= t_min_sim{data_id});
           t_max_calibrated{data_id} = t_domain_max{data_id}(t_domain_max{data_id} <= t_max_sim{data_id});
           
                %IF statements for constraining x-axes
                if numel(t_min_calibrated{data_id}) == 2 
                    
                   t_min_calibrated{data_id} = t_min_shifted{data_id};
                   
                end
                
                if numel(t_max_calibrated{data_id}) == 2
                    
                   t_max_calibrated{data_id} = t_max_shifted{data_id};
                   
                end
                
                if t_min_calibrated{data_id} > t_max_sim{data_id} || t_max_calibrated{data_id} < t_min_sim{data_id}
                    
                   display('Calibrated data is out of domain, please reset');
                   return;
                    
                end
            
           faml_plot_figure_calibration(handles.legend_m_disp_force, data_exp, handles.data_disp_force_calibrated, ptitle, xoutput, youtput, t_min_calibrated{data_id}, t_max_calibrated{data_id}, handles);  
           faml_plot_GUI_calibration(handles.legend_m_disp_force, data_exp, handles.data_disp_force_calibrated, ptitle, xoutput, youtput, t_min_sim{data_id}, t_max_sim{data_id}, handles);
    
        end
        
        %Scale in Y axes
        if get(handles.calibration_checkbox_scale, 'Value')
            
           handles.disp_force{data_id}(:,2) = handles.disp_force{data_id}(:,2)*calibrate_value;
           guidata(hObject, handles);
           handles.data_disp_force_calibrated{data_id} = handles.disp_force{data_id};
           guidata(hObject, handles);
           
           handles.net_calibration{3,data_id}(1,3) = handles.net_calibration{3,data_id}(1,3) * calibrate_value;
           guidata(hObject, handles);
           
           handles.legend_m_disp_force{data_id+1} = strcat('T',num2str(data_id), ...
                                                             '(x =   ', num2str(handles.net_calibration{3,data_id}(1,1)), ...
                                                             ' y =   ', num2str(handles.net_calibration{3,data_id}(1,2)), ...
                                                             ' scale =   ', num2str(handles.net_calibration{3,data_id}(1,3)), ')');
                                                                   
           guidata(hObject, handles);
      
           faml_plot_figure_calibration(handles.legend_m_disp_force, data_exp, handles.data_disp_force_calibrated, ptitle, xoutput, youtput, t_min_sim{data_id}, t_max_sim{data_id}, handles);  
           faml_plot_GUI_calibration(handles.legend_m_disp_force, data_exp, handles.data_disp_force_calibrated, ptitle, xoutput, youtput, t_min_sim{data_id}, t_max_sim{data_id}, handles);
            
        end
        
        %Display warning messages when multiple operations are selected
        if get(handles.calibration_checkbox_shift, 'Value') && ...
           get(handles.calibration_checkbox_scale, 'Value') || ...
           get(handles.calibration_checkbox_shift, 'Value') && ...
           get(handles.calibration_checkbox_shift_time, 'Value') || ...
           get(handles.calibration_checkbox_shift_time, 'Value') && ...
           get(handles.calibration_checkbox_scale, 'Value')
            
           display('Multiple operation selected. Please choose either Scale or Shift');
           return;
           
        end
        
        if get(handles.calibrate_force, 'Value') && ...
           get(handles.calibrate_disp, 'Value') || ...
           get(handles.calibrate_force, 'Value') && ...
           get(handles.calibrate_disp_force, 'Value') || ...
           get(handles.calibrate_disp, 'Value') && ...
           get(handles.calibrate_disp_force, 'Value')
            
           display('Multiple variable selected. Please choose one variable only');
           return;
           
        end
     
    end
    

% --- Resets calibration
% --- Executes on button press in calibration_reset.
function calibration_reset_Callback(hObject, eventdata, handles)
% hObject    handle to calibration_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    cla(handles.axes1);
    cla(handles.axes2);
    cla(handles.axes3);
    
    %Bring number of imported simulation data
    sim_index = handles.sim_index - 1;
   
    %Initialise legend matrices
    for i = 1 : sim_index
        
        %Reset legend
        handles.legend_m_force{i+1} = strcat('T', num2str(i), ' original');
        guidata(hObject, handles);
        
        handles.legend_m_disp{i+1} = strcat('T', num2str(i), ' original');
        guidata(hObject, handles);
        
        handles.legend_m_disp_force{i+1} = strcat('T', num2str(i), ' original');
        guidata(hObject, handles);
        
        %Reset net calibration value
        handles.net_calibration{1,i} = zeros(1,3);
        guidata(hObject, handles);
        handles.net_calibration{1,i}(1,3) = 1;
        guidata(hObject, handles);
        
        handles.net_calibration{2,i} = zeros(1,3);
        guidata(hObject, handles);
        handles.net_calibration{2,i}(1,3) = 1;
        guidata(hObject, handles);
        
        handles.net_calibration{3,i} = zeros(1,3);
        guidata(hObject, handles);
        handles.net_calibration{3,i}(1,3) = 1;
        guidata(hObject, handles);
       
        %Reset data
        handles.data_sim{i} = handles.data_sim_original{i};
        guidata(hObject, handles);
    
        handles.time_force{i} = handles.time_force_original{i};
        guidata(hObject, handles);
    
        handles.time_disp{i} = handles.time_disp_original{i};
        guidata(hObject, handles);
    
        handles.disp_force{i} = handles.disp_force_original{i};
        guidata(hObject, handles);
        
    end
    
    data_exp = handles.data_exp;
    
    %t_min_sim = handles.t_min_sim;
    %t_max_sim = handles.t_max_sim;
    
    %Plot the original experimental and simulation data
    faml_plot_figure_reset(handles.legend_m_force, data_exp, handles.data_sim, '', '');  
    faml_plot_GUI_reset(handles.legend_m_force, data_exp, handles.data_sim, '', '', handles);



%DEFINE GUI COMPONENTS: DO NOT MAKE CHANGES
%==================================================================================%
function calibration_value_Callback(hObject, eventdata, handles)
% hObject    handle to calibration_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of calibration_value as text
%        str2double(get(hObject,'String')) returns contents of calibration_value as a double
    
% --- Executes during object creation, after setting all properties.
function calibration_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calibration_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
    
end

% --- Executes on button press in calibrate_force.
function calibrate_force_Callback(hObject, eventdata, handles)
% hObject    handle to calibrate_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of calibrate_force

% --- Executes on button press in calibrate_disp.
function calibrate_disp_Callback(hObject, eventdata, handles)
% hObject    handle to calibrate_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of calibrate_disp

% --- Executes on button press in calibrate_disp_force.
function calibrate_disp_force_Callback(hObject, eventdata, handles)
% hObject    handle to calibrate_disp_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of calibrate_disp_force

% --- Executes on button press in calibration_checkbox_shift.
function calibration_checkbox_shift_Callback(hObject, eventdata, handles)
% hObject    handle to calibration_checkbox_shift (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of calibration_checkbox_shift

% --- Executes on button press in calibration_checkbox_shift_time.
function calibration_checkbox_shift_time_Callback(hObject, eventdata, handles)
% hObject    handle to calibration_checkbox_shift_time (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of calibration_checkbox_shift_time

% --- Executes on button press in calibration_checkbox_scale.
function calibration_checkbox_scale_Callback(hObject, eventdata, handles)
% hObject    handle to calibration_checkbox_scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of calibration_checkbox_scale

% --- Executes on button press in T1.
function T1_Callback(hObject, eventdata, handles)
% hObject    handle to T1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of T1

% --- Executes on button press in T2.
function T2_Callback(hObject, eventdata, handles)
% hObject    handle to T2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of T2

% --- Executes on button press in T3.
function T3_Callback(hObject, eventdata, handles)
% hObject    handle to T3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of T3

% --- Executes on button press in T4.
function T4_Callback(hObject, eventdata, handles)
% hObject    handle to T4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of T4



