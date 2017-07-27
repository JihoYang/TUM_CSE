
function faml_plot_GUI_calibration(legend_m, data, data_calibrated, ptitle, xoutput, youtput, t_min, t_max, handles)
    
    dim = size(data_calibrated);
    num_sim = dim(2);
    
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
    
    for i = 1 : num_sim
        
        plot(handles.(['axes',num2str(k)]), data_calibrated{i}(:,1), data_calibrated{i}(:,2), 'color', rand(1,3));
        
    end
    
    title(handles.(['axes',num2str(k)]), ptitle(k));
    xlabel(handles.(['axes',num2str(k)]), xoutput(k));
    ylabel(handles.(['axes',num2str(k)]), youtput(k));
    legend(handles.(['axes',num2str(k)]), legend_m, 'Location', 'northeast');
        
    box(handles.(['axes',num2str(k)]), 'on');
    
    xlim(handles.(['axes',num2str(k)]), [t_min, t_max]); 
        
    set(handles.(['axes',num2str(k)]), 'FontSize', 10);
    set(handles.(['axes',num2str(k)]), 'FontName', 'Arial');
    hold(handles.(['axes',num2str(k)]),'on')
    
        if k == 3
            
           cla(handles.(['axes',num2str(k)]));
            
           plot(handles.(['axes',num2str(k)]), data(:,2), data(:,3), 'b'); 
           hold(handles.(['axes',num2str(k)]),'on')
           
           for i = 1:num_sim
               
                plot(handles.(['axes',num2str(k)]), data_calibrated{i}(:,1), data_calibrated{i}(:,2), 'color', rand(1,3));
            
           end
           
           title(handles.(['axes',num2str(k)]), ptitle(k));
           xlabel(handles.(['axes',num2str(k)]), xoutput(k));
           ylabel(handles.(['axes',num2str(k)]), youtput(k));
           legend(handles.(['axes',num2str(k)]), legend_m, 'Location', 'northeast');
           xlim(handles.axes3, 'auto');
           
           box(handles.(['axes',num2str(k)]), 'on');
  
           set(handles.(['axes',num2str(k)]), 'FontSize', 10);
           set(handles.(['axes',num2str(k)]), 'FontName', 'Arial');
           hold(handles.(['axes',num2str(k)]),'on')
           
        end         
        
end