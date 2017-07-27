function faml_plot_GUI_reset(legend_m, data, data_calibrated, t_min, t_max, handles)
        
        dim = size(data_calibrated);
        num_sim = dim(2);
        
        cla(handles.axes1);
           
        plot(handles.axes1, data(:,1), data(:,2), 'b'); 
        hold(handles.axes1,'on')
        
        for i = 1:num_sim
            
            plot(handles.axes1, data_calibrated{i}(:,1), data_calibrated{i}(:,2), 'color', rand(1,3));
            hold on;
        
        end
        
        title(handles.axes1, 'Force vs Time (Exp vs Sim (reset))');
        xlabel(handles.axes1, 'Time (sec)');
        ylabel(handles.axes1, 'Force (kN)');
        legend(handles.axes1, legend_m, 'Location', 'northeast');
        
        box(handles.axes1, 'on');
        %xlim(handles.axes1, [t_min, t_max]); 
        set(handles.axes1, 'FontSize', 10);
        set(handles.axes1, 'FontName', 'Arial');
        
        
        cla(handles.axes2);
           
        plot(handles.axes2, data(:,1), data(:,3), 'b'); 
        hold(handles.axes2,'on')
        
        for i = 1:num_sim
            
            plot(handles.axes2, data_calibrated{i}(:,1), data_calibrated{i}(:,3), 'color', rand(1,3));
            
        end
        
        title(handles.axes2, 'Displacement vs Time (Exp vs Sim (reset))');
        xlabel(handles.axes2, 'Time (sec)');
        ylabel(handles.axes2, 'Displacement (mm)');
        legend(handles.axes2, legend_m, 'Location', 'northeast');
        
        box(handles.axes2, 'on');
        %xlim(handles.axes2, [t_min, t_max]); 
        set(handles.axes2, 'FontSize', 10);
        set(handles.axes2, 'FontName', 'Arial');
        
        cla(handles.axes3);
           
        plot(handles.axes3, data(:,2), data(:,3), 'b'); 
        hold(handles.axes3, 'on')
        
        for i = 1:num_sim
            
            plot(handles.axes3, data_calibrated{i}(:,2), data_calibrated{i}(:,3), 'color', rand(1,3));
            
        end
        
        title(handles.axes3, 'Displacement vs Force (Exp vs Sim (reset))');
        xlabel(handles.axes3, 'Force (kN)');
        ylabel(handles.axes3, 'Diplacement (mm)');
        legend(handles.axes3, legend_m, 'Location', 'northeast');
        
        box(handles.axes3, 'on');
        %xlim(handles.axes3, 'auto'); 
        set(handles.axes3, 'FontSize', 10);
        set(handles.axes3, 'FontName', 'Arial');            
        
end