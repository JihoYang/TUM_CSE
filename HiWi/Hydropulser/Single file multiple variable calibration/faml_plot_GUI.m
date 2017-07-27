
%FAML_PLOT_GUI generates plots for imported data at the pre-defined
%axes in the GUI environment

function faml_plot_GUI(legend_m_sim, data, variables, unit, t_min, t_max, data_type, handles)
    
    
    if data_type == 'Experiment';
        
        for k = 1:2
            
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
     
        plot(handles.axes4, data(:,1), data(:,2));
        title(handles.axes4, strcat(variables(1), ' variation over time'));
        xlabel(handles.axes4, 'Time (sec)');
        ylabel(handles.axes4, strcat(variables(1), unit(1)));
        legend(handles.axes4, legend_m_sim{:,1}, 'Location', 'northeast');
        box(handles.axes4, 'on');
        xlim(handles.axes4, [t_min, t_max]); 
        set(handles.axes4, 'FontSize', 10);
        set(handles.axes4, 'FontName', 'Arial');
        hold(handles.axes4, 'on');
        
        plot(handles.axes5, data(:,1), data(:,3));
        title(handles.axes5, strcat(variables(2), ' variation over time'));
        xlabel(handles.axes5, 'Time (sec)');
        ylabel(handles.axes5, strcat(variables(2), unit(2)));
        legend(handles.axes5, legend_m_sim{:,2}, 'Location', 'northeast');
        box(handles.axes5, 'on');
        xlim(handles.axes5, [t_min, t_max]); 
        set(handles.axes5, 'FontSize', 10);
        set(handles.axes5, 'FontName', 'Arial');
        hold(handles.axes5, 'on');

        plot(handles.axes6, data(:,2), data(:,3));
        
        title(handles.axes6, 'Displacement vs Force');
        xlabel(handles.axes6, 'Force (kN)');
        ylabel(handles.axes6, 'Displacement (mm)');
        legend(handles.axes6, legend_m_sim{:,3}, 'Location', 'northeast');
        
        box(handles.axes6, 'on');
        xlim(handles.axes6, [t_min, t_max]); 
        set(handles.axes6, 'FontSize', 10);
        set(handles.axes6, 'FontName', 'Arial');
        
        hold(handles.axes6, 'on')
        
        plot(handles.axes7, data(:,1), data(:,4));
        box(handles.axes7, 'on');
        legend(handles.axes7, handles.legend_m_sim{:,4}, 'Location', 'northeast');
        
        hold(handles.axes7, 'on')
        
        plot(handles.axes8, data(:,1), data(:,5));
        box(handles.axes8, 'on');
        legend(handles.axes8, handles.legend_m_sim{:,5}, 'Location', 'northeast');
        
        hold(handles.axes8, 'on')

    end
        
end






