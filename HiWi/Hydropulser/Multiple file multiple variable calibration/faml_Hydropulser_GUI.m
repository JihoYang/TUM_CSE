
%===========================================================================%
%|                                                                         |%
%|     Graphical User Interface for pre-processing experimental data       |%
%|                                                                         |%
%|     Model: HYDROPULSER                                                  |%
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
%|     Final update: 14th January 2016                                     |%  
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

% Last Modified by GUIDE v2.5 13-Jan-2016 21:27:01

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
     dimension = size(data_numeric);
     t_min = min(data_numeric(:,1));
     t_max = max(data_numeric(:,1));
     
     output = {'Force', 'Displacement', 'Temperature'};
     unit   = {' (kN)', ' (mm)', ' (\circC)'};
     
     %Create handles for sharing data between callbacks
     handles.data_exp = data_numeric;
     guidata(hObject, handles);
     
     %Plot data in GUI environment
     faml_plot_GUI(data_numeric, dimension, output, unit, t_min, t_max, 'Experiment', handles);
 
     %Plot data in separate figure windows
     faml_plot_figure(data_numeric, dimension, t_min, t_max, output, unit, 'Experiment');
     
     %Display full file path directory in GUI environment
     set(handles.experiment_filepath, 'String', fullpathname);

     
% --- Executes on button press in Browse_SimulationData.
% --- Loads and processes simulation data in text file format (.tab)
function Browse_SimulationData_Callback(hObject, eventdata, handles)
% hObject    handle to Browse_SimulationData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %Bring the selected text file (.txt.) and save its file and path names
    [filename_sim pathname_sim] = uigetfile({'.tab'},'File Selector');
    fullpathname_sim = strcat(pathname_sim, filename_sim);

    %Build full file name and read/import it to MATLAB GUI workspace
    loaddata_sim = fullfile(pathname_sim, filename_sim);
    data_numeric_sim = dlmread(loaddata_sim,'\t',7,0);
    
    %Define parameters for data plotting
    dimension_sim = size(data_numeric_sim);
    t_min_sim = min(data_numeric_sim(:,1));
    t_max_sim = max(data_numeric_sim(:,1));

    output_sim = {'Force', 'Displacement'};
    unit_sim   = {' (kN)', ' (mm)'};
    
    %Save the imported simulation data in MATLAB base workspace
    setappdata(handles.Browse_SimulationData,'data_original_sim', data_numeric_sim);
    
    %Create handles for sharing data
    handles.data_sim = data_numeric_sim;
    guidata(hObject, handles);
    
    handles.t_min_sim=t_min_sim;
    guidata(hObject, handles);
     
    handles.t_max_sim=t_max_sim;
    guidata(hObject, handles);
 
    %Plot data in GUI environment
    faml_plot_GUI(data_numeric_sim, dimension_sim, output_sim, unit_sim, t_min_sim, t_max_sim, 'Simulation', handles);

    %Plot data in separate figure windows
    faml_plot_figure(data_numeric_sim, dimension_sim, t_min_sim, t_max_sim, output_sim, unit_sim, 'Simulation');
    
    %Display full file path directory in GUI environment
    set(handles.simulation_filepath, 'String', fullpathname_sim);

    
% --- Imports simulation data (time & force)
% --- Executes on button press in browse_force.
function browse_force_Callback(hObject, eventdata, handles)
% hObject    handle to browse_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %Bring the selected text file (.txt.) and save its file and path names
     [filename pathname] = uigetfile({'.txt'},'File Selector');
     fullpathname = strcat(pathname, filename);
     
     %Build full file name and read/import it to MATLAB GUI workspace
     loaddata = fullfile(pathname, filename);
     data_numeric_force = dlmread(loaddata,'\t');
     
     %============================================================================%
     %The following lines are used for importing text files with ',' as
     %decimal place indicator
     
     %Reads ',' as strings and converts them into '.'
     
     %data_numeric = [];
 
     %   for i=1:size(data_string,1)
         
     %       data_numeric(i,:) = str2num(char(strrep(data_string(i,:),',','.')'));
          
     %   end
     %============================================================================%
    
     %Create handles for sharing data between callbacks
     handles.data_force_sim = data_numeric_force;
     guidata(hObject, handles);
     
     handles.data_numeric_force = data_numeric_force;
     guidata(hObject, handles);
     
     %Plot data in GUI environment
     cla(handles.axes1);
     plot(handles.axes1, handles.data_exp(:,1), handles.data_exp(:,2), 'b');
     hold on;
     plot(handles.axes1, data_numeric_force(:,1), data_numeric_force(:,2), 'r');
     title(handles.axes1, 'Force vs Time (Exp vs Sim (original))');
     xlabel(handles.axes1, 'Time (sec)');
     ylabel(handles.axes1, 'Force (kN)');
     xlim(handles.axes1, [handles.t_min_sim, handles.t_max_sim]);
     legend(handles.axes1, 'Experiment', 'Simulation','Location','northeast');
     set(handles.axes1, 'FontSize', 10);
     set(handles.axes1, 'FontName', 'Arial');
     
     
     %Plot data in a separate figure
     figure(7)
     h = figure(7);
     cla(h);
     
     plot(handles.data_exp(:,1), handles.data_exp(:,2), 'b');
     hold on;
     
     plot(data_numeric_force(:,1), data_numeric_force(:,2), 'r');
     
     title('Force vs Time (Experiment vs original Simulation)');
     xlabel('Time (sec)');
     ylabel('Force (kN)');
     xlim([handles.t_min_sim, handles.t_max_sim]);
     legend('Experiment', 'Simulation', 'Location', 'northeast');
     box on; grid on;
     set(gca, 'FontSize', 15);
     set(gca, 'FontName', 'Arial');
     
     %Display full file path directory in GUI environment
     set(handles.force_calibration_path, 'String', fullpathname);

     
% --- Imports simulation data (time & displacement)
% --- Executes on button press in browse_disp.
function browse_disp_Callback(hObject, eventdata, handles)
% hObject    handle to browse_disp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Bring the selected text file (.txt.) and save its file and path names

     [filename pathname] = uigetfile({'.txt'},'File Selector');
     fullpathname = strcat(pathname, filename);
     
     %Build full file name and read/import it to MATLAB GUI workspace
     loaddata = fullfile(pathname, filename);
     data_numeric_disp = dlmread(loaddata,'\t');
     
     %============================================================================%
     %The following lines are used for importing text files with ',' as
     %decimal place indicator
     
     %Reads ',' as strings and converts them into '.'
     
     %data_numeric = [];
 
     %   for i=1:size(data_string,1)
         
     %       data_numeric(i,:) = str2num(char(strrep(data_string(i,:),',','.')'));
          
     %   end
     %============================================================================%
    
     %Create handles for sharing data between callbacks
     handles.data_disp_sim = data_numeric_disp;
     guidata(hObject, handles);
     
     handles.data_numeric_disp = data_numeric_disp;
     guidata(hObject, handles);
     
     %Plot data in GUI environment
     cla(handles.axes2);
     plot(handles.axes2, handles.data_exp(:,1), handles.data_exp(:,3), 'b');
     hold on;
     plot(handles.axes2, data_numeric_disp(:,1), data_numeric_disp(:,2), 'r');
     xlim(handles.axes2, [handles.t_min_sim, handles.t_max_sim]);
     title(handles.axes2, 'Disp vs Time (Exp vs Sim (origianl))');
     xlabel(handles.axes2, 'Time (sec)');
     ylabel(handles.axes2, 'Displacement (mm)');
     legend(handles.axes2, 'Experiment', 'Simulation','Location','northeast');
     set(handles.axes2, 'FontSize', 10);
     set(handles.axes2, 'FontName', 'Arial');
     
     %Plot data in a separate figure
     figure(8)
     h = figure(8);
     cla(h);
     
     plot(handles.data_exp(:,1), handles.data_exp(:,3), 'b');
     hold on;
     
     plot(data_numeric_disp(:,1), data_numeric_disp(:,2), 'r');
     
     title('Displacement vs Time (Experiment vs original Simulation)');
     xlabel('Time (sec)');
     ylabel('Displacement (mm)');
     xlim([handles.t_min_sim, handles.t_max_sim]);
     legend('Experiment', 'Simulation', 'Location', 'northeast');
     box on; grid on;
     set(gca, 'FontSize', 15);
     set(gca, 'FontName', 'Arial');
     
     %Display full file path directory in GUI environment
     set(handles.disp_calibration_path, 'String', fullpathname);

% --- Imports simulation data (force & displacement)
% --- Executes on button press in browse_disp_force.
function browse_disp_force_Callback(hObject, eventdata, handles)
% hObject    handle to browse_disp_force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%Bring the selected text file (.txt.) and save its file and path names
     [filename pathname] = uigetfile({'.txt'},'File Selector');
     fullpathname = strcat(pathname, filename);
     
     %Build full file name and read/import it to MATLAB GUI workspace
     loaddata = fullfile(pathname, filename);
     data_numeric_disp_force = dlmread(loaddata,'\t');
     
     handles.data_numeric_disp_force = data_numeric_disp_force;
     
     %============================================================================%
     %The following lines are used for importing text files with ',' as
     %decimal place indicator
     
     %Reads ',' as strings and converts them into '.'
     
     %data_numeric = [];
 
     %   for i=1:size(data_string,1)
         
     %       data_numeric(i,:) = str2num(char(strrep(data_string(i,:),',','.')'));
          
     %   end
     %============================================================================%
    
     %Create handles for sharing data between callbacks
     handles.data_disp_force_sim = data_numeric_disp_force;
     guidata(hObject, handles);
     
     handles.data_numeric_disp_force = data_numeric_disp_force;
     guidata(hObject, handles);
     
     %Plot data in GUI environment
     cla(handles.axes3);
     plot(handles.axes3, handles.data_exp(:,2), handles.data_exp(:,3), 'b');
     hold on;
     plot(handles.axes3, data_numeric_disp_force(:,1), data_numeric_disp_force(:,2), 'r');
     title(handles.axes3, 'Disp vs Force (Exp vs Sim (original))');
     xlabel(handles.axes3, 'Force (kN)');
     ylabel(handles.axes3, 'Displacement (mm)');
     legend(handles.axes3, 'Experiment', 'Simulation','Location','northeast');
     set(handles.axes3, 'FontSize', 10);
     set(handles.axes3, 'FontName', 'Arial');
     
     %Plot data in a separate figure
     figure(9);
     h = figure(9);
     cla(h);
     plot(handles.data_exp(:,2), handles.data_exp(:,3), 'b');
     hold on; 
     plot(data_numeric_disp_force(:,1), data_numeric_disp_force(:,2), 'r');
     title('Displacement vs Force (Experiment vs original Simulation)');
     xlabel('Force (kN)');
     ylabel('Displacement (mm)');
     legend('Experiment', 'Simulation', 'Location', 'northeast');
     box on; grid on;
     set(gca, 'FontSize', 15);
     set(gca, 'FontName', 'Arial');
     
     %Display full file path directory in GUI environment
     set(handles.disp_force_calibration_path, 'String', fullpathname);
     

% --- Calibrates experimental data (displacement only)
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
    
    %Bring computed minimum and maximum time from the simulation dataset
    t_min_sim = handles.t_min_sim;
    t_max_sim = handles.t_max_sim;
    
    %Define strings for data plotting
    ptitle = {'Force vs time (calibrated)', 'Displacement vs time (calibrated)', ...
               'Displacement vs Force (calibrated)'}; 
    xoutput = {'Time (sec)', 'Time (sec)', 'Force (kN)'};
    youtput = {'Force (kN)', 'Displacement (mm)', 'Displacement (mm)'};
    
    %Calibrate force vs time
    if get(handles.calibrate_force, 'Value')
       
       %Condition statements for shift/scale operations
    
        %Shift in Y axes
        if get(handles.calibration_checkbox_shift, 'Value')
           
           handles.data_force_sim(:,2) = handles.data_force_sim(:,2) + calibrate_value;
           guidata(hObject, handles);
           data_force_calibrated = handles.data_force_sim;
            
           faml_plot_figure_calibration(data_exp, data_force_calibrated, ptitle, xoutput, youtput, t_min_sim, t_max_sim, handles);  
           faml_plot_GUI_calibration(data_exp, data_force_calibrated, ptitle, xoutput, youtput, t_min_sim, t_max_sim, handles);
    
        end
        
        %Shift in X axes
        if get(handles.calibration_checkbox_shift_time, 'Value')
           
           handles.data_force_sim(:,1) = handles.data_force_sim(:,1) + calibrate_value;
           guidata(hObject, handles);
           data_force_calibrated = handles.data_force_sim;
           
           t_min_shifted = min(data_force_calibrated(:,1));
           t_max_shifted = max(data_force_calibrated(:,1));
           
           t_domain_min = [t_min_sim, t_min_shifted];
           t_domain_max = [t_max_sim, t_max_shifted];
           
           t_min_calibrated = t_domain_min(t_domain_min >= t_min_sim);
           t_max_calibrated = t_domain_max(t_domain_max <= t_max_sim);
           
                %IF statements for constraining x-axes
                if numel(t_min_calibrated) == 2 
                    
                   t_min_calibrated = t_min_shifted;
                   
                end
                
                if numel(t_max_calibrated) == 2
                    
                   t_max_calibrated = t_max_shifted;
                   
                end
                
                if t_min_calibrated > t_max_sim || t_max_calibrated < t_min_sim
                    
                   display('Calibrated data is out of domain, please reset');
                   return;
                    
                end
            
           faml_plot_figure_calibration(data_exp, data_force_calibrated, ptitle, xoutput, youtput, t_min_calibrated, t_max_calibrated, handles);  
           faml_plot_GUI_calibration(data_exp, data_force_calibrated, ptitle, xoutput, youtput, t_min_sim, t_max_sim, handles);
    
        end
        
        %Scale in Y axes
        if get(handles.calibration_checkbox_scale, 'Value')
            
           handles.data_force_sim(:,2) = handles.data_force_sim(:,2)*calibrate_value;
           guidata(hObject, handles);
           data_force_calibrated = handles.data_force_sim;
      
           faml_plot_figure_calibration(data_exp, data_force_calibrated, ptitle, xoutput, youtput, t_min_sim, t_max_sim, handles);  
           faml_plot_GUI_calibration(data_exp, data_force_calibrated, ptitle, xoutput, youtput, t_min_sim, t_max_sim, handles);
            
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
           
           handles.data_disp_sim(:,2) = handles.data_disp_sim(:,2) + calibrate_value;
           guidata(hObject, handles);
           data_disp_calibrated = handles.data_disp_sim;
            
           faml_plot_figure_calibration(data_exp, data_disp_calibrated, ptitle, xoutput, youtput, t_min_sim, t_max_sim, handles);  
           faml_plot_GUI_calibration(data_exp, data_disp_calibrated, ptitle, xoutput, youtput, t_min_sim, t_max_sim, handles);
    
        end
        
        %Shift in X axes
        if get(handles.calibration_checkbox_shift_time, 'Value')
           
           handles.data_disp_sim(:,1) = handles.data_disp_sim(:,1) + calibrate_value;
           guidata(hObject, handles);
           data_disp_calibrated = handles.data_disp_sim;
           
           t_min_shifted = min(data_disp_calibrated(:,1));
           t_max_shifted = max(data_disp_calibrated(:,1));
           
           t_domain_min = [t_min_sim, t_min_shifted];
           t_domain_max = [t_max_sim, t_max_shifted];
           
           t_min_calibrated = t_domain_min(t_domain_min >= t_min_sim);
           t_max_calibrated = t_domain_max(t_domain_max <= t_max_sim);
           
                %IF statements for constraining x-axes
                if numel(t_min_calibrated) == 2 
                    
                   t_min_calibrated = t_min_shifted;
                   
                end
                
                if numel(t_max_calibrated) == 2
                    
                   t_max_calibrated = t_max_shifted;
                   
                end
                
                if t_min_calibrated > t_max_sim || t_max_calibrated < t_min_sim
                    
                   display('Calibrated data is out of domain, please reset');
                   return;
                    
                end
            
           faml_plot_figure_calibration(data_exp, data_disp_calibrated, ptitle, xoutput, youtput, t_min_calibrated, t_max_calibrated, handles);  
           faml_plot_GUI_calibration(data_exp, data_disp_calibrated, ptitle, xoutput, youtput, t_min_sim, t_max_sim, handles);
    
        end
        
        %Scale in Y axes
        if get(handles.calibration_checkbox_scale, 'Value')
            
           handles.data_disp_sim(:,2) = handles.data_disp_sim(:,2)*calibrate_value;
           guidata(hObject, handles);
           data_disp_calibrated = handles.data_disp_sim;
      
           faml_plot_figure_calibration(data_exp, data_disp_calibrated, ptitle, xoutput, youtput, t_min_sim, t_max_sim, handles);  
           faml_plot_GUI_calibration(data_exp, data_disp_calibrated, ptitle, xoutput, youtput, t_min_sim, t_max_sim, handles);
            
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
           
           handles.data_disp_force_sim(:,2) = handles.data_disp_force_sim(:,2) + calibrate_value;
           guidata(hObject, handles);
           data_disp_force_calibrated = handles.data_disp_force_sim;
            
           faml_plot_figure_calibration(data_exp, data_disp_force_calibrated, ptitle, xoutput, youtput, t_min_sim, t_max_sim, handles);  
           faml_plot_GUI_calibration(data_exp, data_disp_force_calibrated, ptitle, xoutput, youtput, t_min_sim, t_max_sim, handles);
    
        end
        
        %Shift in X axes
        if get(handles.calibration_checkbox_shift_time, 'Value')
           
           handles.data_disp_force_sim(:,1) = handles.data_disp_force_sim(:,1) + calibrate_value;
           guidata(hObject, handles);
           data_disp_force_calibrated = handles.data_disp_force_sim;
           
           t_min_shifted = min(data_disp_force_calibrated(:,1));
           t_max_shifted = max(data_disp_force_calibrated(:,1));
           
           t_domain_min = [t_min_sim, t_min_shifted];
           t_domain_max = [t_max_sim, t_max_shifted];
           
           t_min_calibrated = t_domain_min(t_domain_min >= t_min_sim);
           t_max_calibrated = t_domain_max(t_domain_max <= t_max_sim);
           
                %IF statements for constraining x-axes
                if numel(t_min_calibrated) == 2 
                    
                   t_min_calibrated = t_min_shifted;
                   
                end
                
                if numel(t_max_calibrated) == 2
                    
                   t_max_calibrated = t_max_shifted;
                   
                end
                
                if t_min_calibrated > t_max_sim || t_max_calibrated < t_min_sim
                    
                   display('Calibrated data is out of domain, please reset');
                   return;
                    
                end
            
           faml_plot_figure_calibration(data_exp, data_disp_force_calibrated, ptitle, xoutput, youtput, t_min_calibrated, t_max_calibrated, handles);  
           faml_plot_GUI_calibration(data_exp, data_disp_force_calibrated, ptitle, xoutput, youtput, t_min_sim, t_max_sim, handles);
    
        end
        
        %Scale in Y axes
        if get(handles.calibration_checkbox_scale, 'Value')
            
           handles.data_disp_force_sim(:,2) = handles.data_disp_force_sim(:,2)*calibrate_value;
           guidata(hObject, handles);
           data_disp_force_calibrated = handles.data_disp_force_sim;
      
           faml_plot_figure_calibration(data_exp, data_disp_force_calibrated, ptitle, xoutput, youtput, t_min_sim, t_max_sim, handles);  
           faml_plot_GUI_calibration(data_exp, data_disp_force_calibrated, ptitle, xoutput, youtput, t_min_sim, t_max_sim, handles);
            
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


% --- Resets calibration
% --- Executes on button press in calibration_reset.
function calibration_reset_Callback(hObject, eventdata, handles)
% hObject    handle to calibration_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    

    cla(handles.axes1);
    cla(handles.axes2);
    cla(handles.axes3);
    
    %Bring the calibrated data and reset it as the original imported simulation data
    handles.data_force_sim = handles.data_numeric_force;
    guidata(hObject, handles);
   
    handles.data_disp_sim = handles.data_numeric_disp;
    guidata(hObject, handles);
    
    handles.data_disp_force_sim = handles.data_numeric_disp_force;
    guidata(hObject, handles);
    
    data_exp = handles.data_exp;
    data_sim_original(:,1) = handles.data_force_sim(:,1);
    data_sim_original(:,2) = handles.data_force_sim(:,2);
    data_sim_original(:,3) = handles.data_disp_sim(:,2);
    
    t_min_sim = handles.t_min_sim;
    t_max_sim = handles.t_max_sim;
    
    %Plot the original experimental and simulation data
    faml_plot_figure_reset(data_exp, data_sim_original, t_min_sim, t_max_sim);  
    faml_plot_GUI_reset(data_exp, data_sim_original, t_min_sim, t_max_sim, handles);

%==================================================================================%
