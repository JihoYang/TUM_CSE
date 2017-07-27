
function faml_plot_GUI_calibration(data, data_calibrated, ptitle, xoutput, youtput, t_min, t_max, handles)
            
    if get(handles.calibrate_force, 'Value')
       
       k = 1;
     
    end
    
    if get(handles.calibrate_disp, 'Value')
       
       k = 2;
     
    end
    
    if get(handles.calibrate_disp_force, 'Value')
       
       k = 3;
     
    end
    
    cla(handles.(['axes',num2str(k)]));
            
    plot(handles.(['axes',num2str(k)]), data(:,1), data(:,k+1), 'b'); 
    hold(handles.(['axes',num2str(k)]),'on')
    plot(handles.(['axes',num2str(k)]), data_calibrated(:,1), data_calibrated(:,2), 'r');
        
    title(handles.(['axes',num2str(k)]), ptitle(k));
    xlabel(handles.(['axes',num2str(k)]), xoutput(k));
    ylabel(handles.(['axes',num2str(k)]), youtput(k));
    legend(handles.(['axes',num2str(k)]), 'Experiment', 'Simulation (Calibrated)', 'Location', 'northeast');
        
    box(handles.(['axes',num2str(k)]), 'on');
    
    xlim(handles.(['axes',num2str(k)]), [t_min, t_max]); 
        
    set(handles.(['axes',num2str(k)]), 'FontSize', 10);
    set(handles.(['axes',num2str(k)]), 'FontName', 'Arial');
    
        if k == 3
            
           cla(handles.(['axes',num2str(k)]));
            
           plot(handles.(['axes',num2str(k)]), data(:,2), data(:,3), 'b'); 
           hold(handles.(['axes',num2str(k)]),'on')
           plot(handles.(['axes',num2str(k)]), data_calibrated(:,1), data_calibrated(:,2), 'r');
        
           title(handles.(['axes',num2str(k)]), ptitle(k));
           xlabel(handles.(['axes',num2str(k)]), xoutput(k));
           ylabel(handles.(['axes',num2str(k)]), youtput(k));
           legend(handles.(['axes',num2str(k)]), 'Experiment', 'Simulation (Calibrated)', 'Location', 'northeast');
           xlim(handles.axes3, 'auto');
           
           box(handles.(['axes',num2str(k)]), 'on'); 
            
           set(handles.(['axes',num2str(k)]), 'FontSize', 10);
           set(handles.(['axes',num2str(k)]), 'FontName', 'Arial');
           
        end         
        
end
