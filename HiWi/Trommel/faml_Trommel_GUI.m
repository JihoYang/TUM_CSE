
%===========================================================================%
%|                                                                         |%
%|     Graphical User Interface for post-processing experimental data      |%
%|                                                                         |%
%|     Model: TROMMEL                                                      |%
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
%|     Final update: 02th June 2016                                        |%  
%|                                                                         |%
%===========================================================================%

%TODO:

%1. After remove last row function both the legend and values for each force doesn't show up (only the last force value shows up)


%Builds basic MATLAB GUI environment. DO NOT MAKE CHANGES
%==================================================================================%
function varargout = faml_Trommel_GUI(varargin)
% FAML_Trommel_GUI MATLAB code for faml_Trommel_GUI.fig
%      FAML_Trommel_GUI, by itself, creates a new FAML_Trommel_GUI or raises the existing
%      singleton*.
%
%      H = FAML_Trommel_GUI returns the handle to a new FAML_Trommel_GUI or the handle to
%      the existing singleton*.
%
%      FAML_Trommel_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FAML_Trommel_GUI.M with the given input arguments.
%
%      FAML_Trommel_GUI('Property','Value',...) creates a new FAML_Trommel_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before faml_Trommel_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to faml_Trommel_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help faml_Trommel_GUI

% Last Modified by GUIDE v2.5 02-Jun-2016 18:52:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @faml_Trommel_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @faml_Trommel_GUI_OutputFcn, ...
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


% --- Executes just before faml_Trommel_GUI is made visible.
function faml_Trommel_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to faml_Trommel_GUI (see VARARGIN)

% Choose default command line output for faml_Trommel_GUI
handles.output = hObject;

%Displays FML logo
axes(handles.logo);
imshow('Logo.png');

%Displays photo of the experimental setup 
axes(handles.photo);
imshow('Photo.png');

% Update handles structure
guidata(hObject, handles);
% UIWAIT makes faml_Trommel_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%Create a matrix for saving velocity and damping coefficients
damping_m= zeros(1,4);
handles.damping_m = damping_m;
guidata(hObject, handles);

set(handles.veldamp,'data',[])


% --- Outputs from this function are returned to the command line.
function varargout = faml_Trommel_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
%==================================================================================%



%FOLLOWING FUNCTIONS MAY BE MODIFIED
%==================================================================================%

% --- Executes on button press in Browse_ExperimentalData.
function Browse_ExperimentalData_Callback(hObject, eventdata, handles)
% hObject    handle to Browse_ExperimentalData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %Bring the selected text file (.xlsx) and save its file and path names
    %[filename_exp pathname_exp] = uigetfile({'.xlsx'},'File Selector');
    [filename_exp pathname_exp] = uigetfile({'.asc'},'File Selector');
    fullpathname_exp = strcat(pathname_exp, filename_exp);

    %Build full file name and read/import it to MATLAB GUI workspace
    loaddata_exp = fullfile(pathname_exp, filename_exp);
    data_cellarray = textread(loaddata_exp,'%s','delimiter','\n');
    %data_cellarray = importdata(loaddata_exp);
    
    %Remove header (assuming the first row of data starts from the 39th)
    data_cellarray(1:38) = [];
    
    %Change commas to dots and create numeric matrix
    data_numeric_exp = double(str2num(char(strrep(data_cellarray,',','.')')));
    
    
    %============================================================================%
    %The following lines are used when importing .xlsx files
    
    %Remove unnecessary cells
    %NaN = find(isnan(data_cellarray.data(:,1)));
    %data_cellarray.data(NaN,:) = [];
    %dim_data = size(data_cellarray.data);
    
    %    if dim_data(1) > 10;
            
    %        data_cellarray.data(:,(12:end))= [];
    %        data_cellarray.data(:,1) = [];
            
    %    else
            
    %        return
            
    %    end
    
    %data_numeric_exp = data_cellarray.data;
   %============================================================================%
   
    
    %Define/Compute parameters for data plotting
    freq = 100;
    dimension_exp = size(data_numeric_exp);
    t_exp = zeros(dimension_exp(1),1);
    t_exp(2:end) = 1/freq*(1:dimension_exp(1)-1);           
     
    %Save the imported experimental data in MATLAB base workspace
    setappdata(handles.Browse_ExperimentalData,'data_original_exp', data_numeric_exp);
    
    %Create handles for sharing data
    handles.data_exp = data_numeric_exp;
    guidata(hObject, handles);
    
    handles.t_exp = t_exp;
    guidata(hObject, handles);
    
    handles.freq = freq;
    guidata(hObject, handles);
    
    
    %Plot data in GUI environment
    cla(handles.axes1); cla(handles.axes2); cla(handles.axes3);
    faml_plot_GUI(t_exp, data_numeric_exp, t_exp(1), t_exp(end), 'Experiment', handles);
    
    %Plot data in separate figure windows
    cla(figure(1)); cla(figure(2)); cla(figure(3));
    faml_plot_figure(t_exp, data_numeric_exp, '', t_exp(1), t_exp(end), 'Experiment');
  

%Display full file path directory in GUI environment
set(handles.exp_filepath, 'String', fullpathname_exp);


% --- Executes on button press in Browse_SimulationData.
function Browse_SimulationData_Callback(hObject, eventdata, handles)
% hObject    handle to Browse_SimulationData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

    %Bring the selected text file (.dat) and save its file and path names
    [filename_sim pathname_sim] = uigetfile({'.dat'},'File Selector');
    fullpathname_sim = strcat(pathname_sim, filename_sim);

    %Build full file name and read/import it to MATLAB GUI workspace
    loaddata_sim = fullfile(pathname_sim, filename_sim);
    data_cellarray = importdata(loaddata_sim);
    data_numeric_sim = data_cellarray.data;
    
    %Define matrices & vectors needed for calibration
    t_sim = data_numeric_sim(:,1);
    
    %Create a matrix for saving calibration values
    net_calibration = zeros(1,3);
    net_calibration(1,3) = 1;
    
    %Create a matrix for showing the net calibration value
    legend_m_defl = cell(2,1);
    legend_m_defl(1,1) = {'Experiment'};
    
    %Save the imported simulation data in MATLAB base workspace
    setappdata(handles.Browse_SimulationData,'data_original_sim', data_numeric_sim);
     
    %Create handles for sharing data
    handles.data_sim = data_numeric_sim;
    guidata(hObject, handles);
    
    handles.t_sim= t_sim;
    guidata(hObject, handles);
    
    handles.t_min_sim=t_sim(1);
    guidata(hObject, handles);
      
    handles.t_max_sim=t_sim(end);
    guidata(hObject, handles);
    
    handles.net_calibration = net_calibration;
    guidata(hObject, handles);
    
    handles.legend_m_defl = legend_m_defl;
    guidata(hObject, handles);
    
    %Bring the deflection data (experimental)
    defl_exp = handles.defl_exp;
    
    %Define time range for plotting
        if handles.t_exp(1) < t_sim(1);
            
            t_min = t_sim(1);
            
        else
            
            t_min = handles.t_exp(1);
            
        end
        
        if handles.t_exp(end) > t_sim(end);
            
            t_max = t_sim(end);
            
        else
            
            t_max = handles.t_exp(1);
            
        end
        
    %Create handle for t_min,max
    handles.t_min = t_min;
    guidata(hObject, handles);
    
    handles.t_max = t_max;
    guidata(hObject, handles);
  
    %Plot data in GUI environment
    faml_plot_GUI(t_sim, data_numeric_sim, t_min, t_max, ' ', handles);
 
    %Plot data in separate figure windows
    faml_plot_figure(t_sim, data_numeric_sim, defl_exp, t_min, t_max, ' ');
     
    %Display full file path directory in GUI environment
    set(handles.sim_filepath, 'String', fullpathname_sim);


% --- Executes on button press in calibration.
function calibration_Callback(hObject, eventdata, handles)
% hObject    handle to calibration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Bring user specified offset and calibration factor
offset_value = str2double(get(handles.offset, 'String'));
correction_factor = str2double(get(handles.correction_factor, 'String'));

%Bring imported experimental data and time vector
data_exp = handles.data_exp;
t_exp = handles.t_exp;

%Calibrate to calculate the deflection
data_exp(:,4) = (data_exp(:,4)+offset_value).*correction_factor;

%Save the calibrated data (deflection)
defl_exp = [t_exp data_exp(:,4)];
handles.defl_exp = defl_exp;
guidata(hObject, handles);

%Plot deflection
faml_plot_GUI_calibration(t_exp, '', defl_exp, t_exp(1), t_exp(end), 'Experiment', '', handles);
faml_plot_figure_calibration(t_exp, '', defl_exp, t_exp(1), t_exp(end), 'Experiment', '', handles);  


% --- Executes on button press in damping.
function calculate_damping_Callback(hObject, eventdata, handles)
% hObject    handle to calculate_damping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Bring the loaded experimental data and define deflection matrix
defl_exp = handles.defl_exp;
t_exp = handles.t_exp;
freq = handles.freq;
vel = str2double(get(handles.velocity, 'String'));
force = str2double(get(handles.force, 'String'));

%Detect peaks
pks = findpeaks(defl_exp(:,2));
number_peak = size(pks);
avg_peak = mean(pks);

    while number_peak(1) > 10
        
        pks = findpeaks(pks);
        number_peak = size(pks);
        avg_peak = mean(pks);
        
    end
    
pks_max = max(pks);
pks_min = min(pks);
pks_diff = abs(pks_max - pks_min);

remove = find(pks < pks_max - 0.9*pks_diff);
pks(remove) = [];

%Remove initial oscillation before first peak
first_peak = find(abs(defl_exp(:,2) - pks(1)) < 1000*eps);
defl_exp((1:first_peak-1),:) = [];

%Compute excitation period
first_peak = find(abs(defl_exp(:,2) - pks(1)) < 1000*eps);
second_peak = find(abs(defl_exp(:,2) - pks(2)) < 1000*eps);
period = t_exp(second_peak) - t_exp(first_peak);
number_excitation = floor((defl_exp(end,1)-defl_exp(1,1))/period);

%Create a matrix for saving each excitation cycle
cycle = find(abs(defl_exp(:,1) - (defl_exp(1,1) + period)) < 1000*eps);

    for i = 1 : cycle - 1
        
        for j = 1 : number_excitation
            
            excit(i,j) = defl_exp((j-1)*(cycle-1)+i,2);
            
        end
        
    end
    
%Calculate the average excitation cycle
excit_tp = excit';
avg_defl_y = mean(excit_tp)';
avg_defl_x = defl_exp((1:cycle-1),1);
avg_defl(:,1) = avg_defl_x;
avg_defl(:,2) = avg_defl_y;

%Remove last bit of uncessary data
avg_defl((end-20:end),:) = [];

%Shift the averaged deflection's time so it starts from t = 0s
min_time_avg_defl = min(avg_defl(:,1));
avg_defl(:,1) = avg_defl(:,1) - min_time_avg_defl;

%Unit conversion from mm to m
%avg_defl(:,2) = avg_defl(:,2)*0.001; 

%Calculate damping coefficient
peak_max = max(avg_defl(:,2));
locs_max = find(peak_max - avg_defl(:,2) < 1000*eps);
[peak_avg locs_avg] = findpeaks(avg_defl(:,2));
dim_peak_avg = size(peak_avg);
peak = zeros(dim_peak_avg(1)+1,1);
locs = zeros(dim_peak_avg(1)+1,1);
peak(1) = peak_max;
peak(2:end) = peak_avg;
locs(1) = locs_max;
locs(2:end) = locs_avg;
P = polyfit(avg_defl(locs(1:5),1), log(peak(1:5)),1); %Only the first 5 peaks are used
delta = -P(1);
mass = force/9.81;
damping_coefficient = 2*delta*mass;

%Damping ratio calculation - OPTIONAL
%dim_avg_defl = size(avg_defl);
%dim_peak = size(peak);
%cycle_damping = dim_peak(1);
%zeta = (log(peak(1)/peak(end)))/((cycle_damping-1)*2*pi);

%Plot deflection with exponential curve


%OPTION 1 (Plot the exponential curve based on 5 first peaks with all the data)
cla(handles.axes5);
plot(handles.axes5, avg_defl(locs,1), peak, 'ro');                                    
title(handles.axes5, 'Exponentially decaying vibration');
xlabel(handles.axes5, 'Time (s)');
ylabel(handles.axes5, 'Deflection (mm)');
box(handles.axes5, 'on');
hold(handles.axes5, 'on');
plot(handles.axes5, avg_defl(:,1), avg_defl(:,2));                                        
plot(handles.axes5, avg_defl(:,1), exp(P(2)+P(1).*avg_defl(:,1)), 'color', [0 0.5 0]);
legend(handles.axes5, 'Peak', 'Averaged deflection', 'Exponential curve', ... 
                      'Location', 'northeast');  
                  
set(handles.axes5, 'FontSize', 10);
set(handles.axes5, 'FontName', 'Arial');

cla(figure(5));
figure(5);
plot(avg_defl(locs,1), peak, 'ro');
hold on;
plot(avg_defl(:,1), avg_defl(:,2));
plot(avg_defl(:,1), exp(P(2)+P(1).*avg_defl(:,1)), 'color', [0 0.5 0]);
title('Exponentially decaying vibration');
xlabel('Time (s)');
ylabel('Deflection (mm)');
box on;
legend('Peak', 'Averaged deflection', 'Exponential curve', ...
        'Location', 'northeast');
set(gca, 'FontName', 'Arial');
set(gca, 'FontSize', 15);

%OPTION 2 (Plot the exponential curve based on 5 first peaks with just first 5 peaks)
%OPTION 2 ALSO CREATES AN ADDITIONAL FIGURE FOR AVERAGED DEFLECTION
% cla(handles.axes5);
% plot(handles.axes5, avg_defl(locs(1:5),1), peak(1:5), 'ro');                          
% title(handles.axes5, 'Exponentially decaying vibration');
% xlabel(handles.axes5, 'Time (s)');
% ylabel(handles.axes5, 'Deflection (mm)');
% box(handles.axes5, 'on');
% hold(handles.axes5, 'on');
% plot(handles.axes5, avg_defl((1:locs(5)),1) , avg_defl((1:locs(5)),2));                 
% plot(handles.axes5, avg_defl((1:locs(5)),1), exp(P(2)+P(1).*avg_defl((1:locs(5)),1))); 
% legend(handles.axes5, 'Peak', 'Averaged deflection', 'Exponential curve', ... 
%                      'Location', 'northeast');  
%                   
% set(handles.axes5, 'FontSize', 10);
% set(handles.axes5, 'FontName', 'Arial');
% cla(figure(5));
% figure(5);
% plot(avg_defl(locs(1:5),1), peak(1:5), 'ro');
% hold on;
% plot(avg_defl((1:locs(5)),1) , avg_defl((1:locs(5)),2));
% plot(avg_defl((1:locs(5)),1), exp(P(2)+P(1).*avg_defl((1:locs(5)),1)));
% title('Exponentially decaying vibration');
% xlabel('Time (s)');
% ylabel('Deflection (mm)');
% box on;
% legend('Peak', 'Averaged deflection', 'Exponential curve', ...
%        'Location', 'northeast');
% set(gca, 'FontName', 'Arial');
% set(gca, 'FontSize', 15);
% cla(figure(6));
% figure(6);
% plot(avg_defl(:,1), avg_defl(:,2));
% title('Averaged Deflection');
% xlabel('Time (s)');
% ylabel('Deflection (mm)');
% box on;
% legend('Averaged Deflection', 'Location', 'northeast');
% set(gca, 'FontName', 'Arial');
% set(gca, 'FontSize', 15);


%Save the computed damping coefficient
handles.damping_m(end+1,:) = [delta, vel, damping_coefficient, force];
guidata(hObject, handles);

%Display the computed data in the GUI environment (table)
set(handles.veldamp, 'data', handles.damping_m);


% --- Executes on button press in plot_damping.
function plot_damping_Callback(hObject, eventdata, handles)
% hObject    handle to plot_damping (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Bring the computed data
damping_m = handles.damping_m;

    if damping_m(1,:) < eps;
        
        %Remove the first row with 0 elements
        damping_m(1,:) = [];
        
        %Save the updated data matrix
        handles.damping_m = damping_m;
        guidata(hObject, handles);
        
        %Find the last row used for the first plot
        dim_damping_m = size(damping_m);
        handles.last_row = dim_damping_m(1);
        guidata(hObject, handles);
        
        %Construct exponential curve for data fitting
        P = polyfit(damping_m(:,2), log(damping_m(:,3)), 1);
        
        %Create an x axis array for plotting the exponential curve
        cf_x = (damping_m(1,2):0.1:damping_m(end,2));
        
        %Extract force value for assigning legends
        legend_f = damping_m(handles.last_row, 4);
        
        %Create and save a matrix for saving the extracted force values
        legend_string = {'Wheel load'};
        legend_m(1) = strcat(legend_string, ' = ', '', num2str(legend_f), 'N');
        legend_m(end+1) = strcat('Exponential curve (', legend_string, ' = ', num2str(legend_f), 'N)');
        handles.legend_m = legend_m;
        guidata(hObject, handles);
        
        %Plot velocity vs damping coefficient in the GUI environment
        cla(handles.axes6);
        plot(handles.axes6, damping_m(:,2), damping_m(:,3), 'o');
        title(handles.axes6, 'Damping coefficient vs Velocity');
        xlabel(handles.axes6, 'Velocity (mm/s)');
        ylabel(handles.axes6, 'Damping coefficient (N\cdots/m)');
       
        box(handles.axes6, 'on');
        set(handles.axes6, 'FontSize', 10);
        set(handles.axes6, 'FontName', 'Arial');
        hold(handles.axes6, 'on');
        
        plot(handles.axes6, cf_x, exp(P(2)+P(1).*cf_x));
        legend(handles.axes6, legend_m(1:end), 'Location', 'northeast');

        %Plot velocity vs damping coefficient in additional figure
        figure(7);
        plot(damping_m(:,2),damping_m(:,3), 'o');
        title('Damping coefficient vs Velocity');
        xlabel('Velocity (mm/s)');
        ylabel('Damping coefficient (N\cdots/m)');
        
        set(gca, 'FontSize', 15);
        set(gca, 'FontName', 'Arial');
        box on;
        hold on;
        
        plot(cf_x, exp(P(2)+P(1).*cf_x));
        
        legend(legend_m(1:end), 'Location', 'northeast');
        
        %Display the updated table in the GUI environment
        set(handles.veldamp, 'data', damping_m);
        
    else
        
        %Bring the necessary data for plotting from the previous IF statement
        last_row = handles.last_row;
        legend_m = handles.legend_m;
        
        %Extract force value for assigning legends
        legend_f = damping_m(last_row+1, 4);
        
        %Construct exponential curve for data fitting
        P = polyfit(damping_m((last_row+1):end,2), log(damping_m((last_row+1:end),3)), 1);
        
        %Create an x axis array for plotting the exponential curve
        cf_x = (damping_m(last_row+1,2):0.1:damping_m(end,2));
        
        %Update and save the matrix for saving the extracted force values
        legend_string = {'Wheel load'};
        legend_m(end+1) = strcat(legend_string, ' = ', '', num2str(legend_f), 'N');
        legend_m(end+1) = strcat('Exponential curve (', legend_string, ' = ', num2str(legend_f), 'N)');
        handles.legend_m = legend_m;
        guidata(hObject, handles);
        
        %Plot velocity vs damping coefficient in the GUI environment
        plot(handles.axes6, damping_m((last_row+1:end),2), damping_m((last_row+1:end),3), 'o');
        
        hold(handles.axes6, 'on');
        
        plot(handles.axes6, cf_x, exp(P(2)+P(1).*cf_x));
        legend(handles.axes6, legend_m, 'Location', 'northeast');
        
        %Plot velocity vs damping coefficient in additional figure
        figure(7);
        plot(damping_m((last_row+1:end),2),damping_m((last_row+1:end),3), 'o');
        
        hold on;
        
        plot(cf_x, exp(P(2)+P(1).*cf_x));
        legend(legend_m, 'Location', 'northeast');
        
        %Update the last row used for the plot
        dim_damping_m = size(damping_m);
        handles.last_row = dim_damping_m(1);
        guidata(hObject, handles);
          
    end

    
% --- Executes on button press in removerow.
function removerow_Callback(hObject, eventdata, handles)
% hObject    handle to removerow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in reset.

    if handles.damping_m(1,:) < eps;

       handles.damping_m(end,:) = [];
       guidata(hObject, handles);

    else

       handles.damping_m(end,:) = [];
       guidata(hObject, handles);

       handles.last_row = handles.last_row - 1;
       guidata(hObject, handles);

   end
   
set(handles.veldamp, 'data', handles.damping_m);     
        
  


function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Reset the matrix for saving velocity and damping coefficients
damping_m= zeros(1,4);
handles.damping_m = damping_m;
guidata(hObject, handles);

%Clear the damping coefficient vs velocity plot
cla(handles.axes5);
cla(handles.axes6);
h1 = figure(6);
h2 = figure(5);
cla(h1);
cla(h2);

set(handles.veldamp,'data',[])


% --- Executes on button press in export_data.
function export_data_Callback(hObject, eventdata, handles)
% hObject    handle to export_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Create a directory for data saving
    if ~exist('Exported Data_MATLAB', 'dir');
    
        mkdir('Exported Data_MATLAB');
    
    end
    
%Bring the computed data
M = handles.damping_m;

%Find the force values
F = unique(M(:,4));

%Find the size of the table
dim_M = size(M);

%Find the index before the force value changes
changedIndex = diff(M(:,4)) ~= 0;

%Find the index from and until which a constant force value is used
F_changed = find(abs(changedIndex(:)-1)<eps);
F_changed(end+1) = dim_M(1);

%Define values for iteration
dim_F = size(F);
last_row = 0;

%Create and export the data matrices
    for i = 1:dim_F(1)

        mat2export{i} = M((last_row+1:F_changed(i)),(1:3));
    
        %Define the filename
        filename = strcat('TP_F_', num2str(F(i)), 'N.txt');
    
        %Export the data matrix
        dlmwrite([pwd '/Exported Data_MATLAB/', filename], mat2export{i}, 'delimiter', '\t');
    
        dim_mat2export = size(mat2export{1,i});
        last_row = dim_mat2export(1) + last_row;
    
    end

    
% --- Executes on button press in calibrate_deflection.
function calibrate_deflection_Callback(hObject, eventdata, handles)
% hObject    handle to calibrate_deflection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Bring user specified calibration value from calibration_value_deflection callback
calibrate_value_defl = str2double(get(handles.calibration_value_deflection, 'String'));

%Bring calibrated deflection data (experimental)
defl_exp = handles.defl_exp;
 
%Bring computed minimum and maximum time 
t_min = handles.t_min;
t_max = handles.t_max;

t_min_sim = handles.t_min_sim;
t_max_sim = handles.t_max_sim;

%Condition statements for shift/scale operations
    
        %Shift in Y axes
        if get(handles.shift_y, 'Value')
           
           handles.data_sim(:,2) = handles.data_sim(:,2) + calibrate_value_defl;
           guidata(hObject, handles);
           data_calibrated = [handles.data_sim(:,1), handles.data_sim(:,2)];
                               
           handles.net_calibration(1,2) = handles.net_calibration(1,2) + calibrate_value_defl;
           guidata(hObject, handles);
        
           handles.legend_m_defl{2,1} = strcat('Simulation (x =   ', num2str(handles.net_calibration(1,1)), ...
                                       ' y =   ', num2str(handles.net_calibration(1,2)), ...
                                       ' scale =   ', num2str(handles.net_calibration(1,3)), ')');
           guidata(hObject, handles);
           
           faml_plot_figure_calibration('', defl_exp, data_calibrated, t_min, t_max, 'Comparison', handles.legend_m_defl, handles);  
           faml_plot_GUI_calibration('', defl_exp, data_calibrated, t_min, t_max, 'Comparison', handles.legend_m_defl, handles);
    
        end
        
        %Shift in X axes
        if get(handles.shift_x, 'Value')
           
           handles.data_sim(:,1) = handles.data_sim(:,1) + calibrate_value_defl;
           guidata(hObject, handles);
           data_calibrated = [handles.data_sim(:,1), handles.data_sim(:,2)];
           
           handles.net_calibration(1,1) = handles.net_calibration(1,1) + calibrate_value_defl;
           guidata(hObject, handles);
        
           handles.legend_m_defl{2,1} = strcat('Simulation (x =   ', num2str(handles.net_calibration(1,1)), ...
                                       ' y =   ', num2str(handles.net_calibration(1,2)), ...
                                       ' scale =   ', num2str(handles.net_calibration(1,3)), ')');
           guidata(hObject, handles);
           
           t_min_shifted = min(data_calibrated(:,1));
           t_max_shifted = max(data_calibrated(:,1));
           
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
            
           faml_plot_figure_calibration('', defl_exp, data_calibrated, t_min_calibrated, t_max_calibrated, 'Comparison', handles.legend_m_defl, handles);  
           faml_plot_GUI_calibration('', defl_exp, data_calibrated, t_min_calibrated, t_max_calibrated, 'Comparison', handles.legend_m_defl, handles);
           
        end
        
        %Scale in Y axes
        if get(handles.scale, 'Value')
            
           handles.data_sim(:,2) = handles.data_sim(:,2)*calibrate_value_defl;
           guidata(hObject, handles);
           data_calibrated = [handles.data_sim(:,1), handles.data_sim(:,2)];
           
           handles.net_calibration(1,3) = handles.net_calibration(1,3)*calibrate_value_defl;
           guidata(hObject, handles);
        
           handles.legend_m_defl{2,1} = strcat('Simulation (x =   ', num2str(handles.net_calibration(1,1)), ...
                                                ' y =   ', num2str(handles.net_calibration(1,2)), ...
                                                ' scale =   ', num2str(handles.net_calibration(1,3)), ')');
           guidata(hObject, handles);
      
           faml_plot_figure_calibration('', defl_exp, data_calibrated, t_min, t_max, 'Comparison', handles.legend_m_defl, handles);  
           faml_plot_GUI_calibration('', defl_exp, data_calibrated, t_min, t_max, 'Comparison', handles.legend_m_defl, handles);
               
           
           
        end
        
        %Display warning messages when multiple operations are selected
        if get(handles.shift_y, 'Value') && ...
           get(handles.scale, 'Value') || ...
           get(handles.shift_y, 'Value') && ...
           get(handles.shift_x, 'Value') || ...
           get(handles.shift_x, 'Value') && ...
           get(handles.scale, 'Value')
            
           display('Multiple operation selected. Please choose either Scale or Shift');
           return;
           
        end

        
% --- Executes on button press in reset_defl.
function reset_defl_Callback(hObject, eventdata, handles)
% hObject    handle to reset_defl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
    
    %Reset the matrices for saving net calibration value and legend
    handles.net_calibration = zeros(1,3);
    handles.net_calibration(1,3) = 1;
    guidata(hObject, handles);
    
    handles.legend_m_defl = cell(2,1);
    handles.legend_m_defl = {'Experiment'};
    guidata(hObject, handles);
    
    %Bring the calibrated data and reset it as the original imported simulation data
    data_original_sim = getappdata(handles.Browse_SimulationData, 'data_original_sim');
    
    handles.data_sim = data_original_sim;
    guidata(hObject, handles);
    
    defl_exp = handles.defl_exp;
    
    t_min = handles.t_min;
    t_max = handles.t_max;
    
    %Plot the original experimental and simulation data (both deflection)
    faml_plot_figure_reset(defl_exp, handles.data_sim, t_min, t_max);  
    faml_plot_GUI_reset(defl_exp, handles.data_sim, t_min, t_max, handles);


%DEFINE GUI COMPONENTS: DO NOT MAKE CHANGES
%==================================================================================%

function exp_filepath_Callback(hObject, eventdata, handles)
% hObject    handle to exp_filepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of exp_filepath as text
%        str2double(get(hObject,'String')) returns contents of exp_filepath as a double

% --- Executes during object creation, after setting all properties.
function exp_filepath_CreateFcn(hObject, eventdata, handles)
% hObject    handle to exp_filepath (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function offset_Callback(hObject, eventdata, handles)
% hObject    handle to offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of offset as text
%        str2double(get(hObject,'String')) returns contents of offset as a double

% --- Executes during object creation, after setting all properties.
function offset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to offset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function correction_factor_Callback(hObject, eventdata, handles)
% hObject    handle to correction_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of correction_factor as text
%        str2double(get(hObject,'String')) returns contents of correction_factor as a double

% --- Executes during object creation, after setting all properties.
function correction_factor_CreateFcn(hObject, eventdata, handles)
% hObject    handle to correction_factor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function velocity_Callback(hObject, eventdata, handles)
% hObject    handle to velocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of velocity as text
%        str2double(get(hObject,'String')) returns contents of velocity as a double

% --- Executes during object creation, after setting all properties.
function velocity_CreateFcn(hObject, eventdata, handles)
% hObject    handle to velocity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function force_Callback(hObject, eventdata, handles)
% hObject    handle to force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of force as text
%        str2double(get(hObject,'String')) returns contents of force as a double

% --- Executes during object creation, after setting all properties.
function force_CreateFcn(hObject, eventdata, handles)
% hObject    handle to force (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in shift_x.
function shift_x_Callback(hObject, eventdata, handles)
% hObject    handle to shift_x (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of shift_x

% --- Executes on button press in shift_y.
function shift_y_Callback(hObject, eventdata, handles)
% hObject    handle to shift_y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of shift_y

% --- Executes on button press in scale.
function scale_Callback(hObject, eventdata, handles)
% hObject    handle to scale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of scale

function calibration_value_deflection_Callback(hObject, eventdata, handles)
% hObject    handle to calibration_value_deflection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of calibration_value_deflection as text
%        str2double(get(hObject,'String')) returns contents of calibration_value_deflection as a double

% --- Executes during object creation, after setting all properties.
function calibration_value_deflection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calibration_value_deflection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%==================================================================================%

