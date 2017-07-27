
%FAML_PLOT_GUI generates plots for imported data at the pre-defined
%axes in the GUI environment

function faml_plot_GUI(data, dimension, variables, unit, t_min, t_max, data_type, handles)
    
    
    if data_type == 'Experiment';
        
        for k = 1:dimension(2)-2
            
            plot(handles.(['axes',num2str(k)]), data(:,1), data(:,k+1));
        
            title(handles.(['axes',num2str(k)]), strcat(variables(k), ' variation over time'));
            xlabel(handles.(['axes',num2str(k)]), 'Time (sec)');
            ylabel(handles.(['axes',num2str(k)]), strcat(variables(k), unit(k)));
            legend(handles.(['axes',num2str(k)]), strcat(variables(k), ' vs Time'), 'Location', 'northeast');
        
            box(handles.(['axes',num2str(k)]), 'on');
            xlim(handles.(['axes',num2str(k)]), [t_min, t_max]); 
            set(handles.(['axes',num2str(k)]), 'FontSize', 10);
            set(handles.(['axes',num2str(k)]), 'FontName', 'Arial');
            
            hold(handles.(['axes',num2str(k)]), 'on')
        
        end
        
        plot(handles.axes3, data(:,2), data(:,3));
        
        title(handles.axes3, 'Displacement variation over time');
        xlabel(handles.axes3, 'Force (kN)');
        ylabel(handles.axes3, 'Displacement (mm)');
        legend(handles.axes3, 'Displacement vs Force');
        
        box(handles.axes3, 'on');
        set(handles.axes3, 'FontSize', 10);
        set(handles.axes3, 'FontName', 'Arial');
        
        hold(handles.axes3, 'on')
        
    else
        
        for k = 1:2
            
            plot(handles.(['axes',num2str(k+3)]), data(:,1), data(:,k+1));
        
            title(handles.(['axes',num2str(k+3)]), strcat(variables(k), ' variation over time'));
            xlabel(handles.(['axes',num2str(k+3)]), 'Time (sec)');
            ylabel(handles.(['axes',num2str(k+3)]), strcat(variables(k), unit(k)));
            legend(handles.(['axes',num2str(k+3)]), strcat(variables(k), ' vs Time'), 'Location', 'northeast');
        
            box(handles.(['axes',num2str(k+3)]), 'on');
            xlim(handles.(['axes',num2str(k+3)]), [t_min, t_max]); 
            set(handles.(['axes',num2str(k+3)]), 'FontSize', 10);
            set(handles.(['axes',num2str(k+3)]), 'FontName', 'Arial');
            
            hold(handles.(['axes',num2str(k+3)]), 'on')
        
        end
        
        plot(handles.axes6, data(:,2), data(:,3));
        
        title(handles.axes6, 'Displacement vs Force');
        xlabel(handles.axes6, 'Force (kN)');
        ylabel(handles.axes6, 'Displacement (mm)');
        legend(handles.axes6, 'Displacement vs Force', 'Location', 'northeast');
        
        box(handles.axes6, 'on'); 
        xlim(handles.axes6, [t_min, t_max]); 
        set(handles.axes6, 'FontSize', 10);
        set(handles.axes6, 'FontName', 'Arial');
        
        hold(handles.axes6, 'on')
        
    end
        
end






