
     %Bring the selected text file (.txt.) and save its file and path names
     [filename pathname] = uigetfile({'.txt'},'File Selector');
     fullpathname = strcat(pathname, filename);
     
     %Build full file name and read/import it to MATLAB GUI workspace
     loaddata = fullfile(pathname, filename);
     data_numeric = dlmread(loaddata,'\t');
     

     %Define parameters for data plotting
     t_min = min(data_numeric(:,1));
     t_max = max(data_numeric(:,1));
     
     output = {'Force', 'Displacement', 'Temperature'};
     unit   = {' (kN)', ' (mm)', ' (\circC)'};
     
     
     for sim_index = 1 : 4
    
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
   
    
    %Initialise calibrated simulation data cell array for calibration
    data_force_calibrated{sim_index} = [time_force{sim_index}(:), data_numeric_sim{sim_index}(:,2)];
    data_disp_calibrated{sim_index} = [time_disp{sim_index}(:), data_numeric_sim{sim_index}(:,2)];
    data_disp_force_calibrated{sim_index} = disp_force{sim_index};
    
    
    %Bring user specified calibration value from calibration_value callback
    calibrate_value = 100;
    
    

    
    
     end
     
    
    
    
    
    
    
    
    
    
    
    
    
    
