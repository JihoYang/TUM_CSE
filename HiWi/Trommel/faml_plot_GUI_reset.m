function faml_plot_GUI_reset(data, data_calibrated, t_min, t_max, handles)
                             
        cla(handles.axes4);
           
        plot(handles.axes4, data(:,1), data(:,2), 'b'); 
        hold(handles.axes4,'on')
        plot(handles.axes4, data_calibrated(:,1), data_calibrated(:,2), 'r');
        
        title(handles.axes4, 'Deflection Comparison (Exp vs Sim (reset))');
        xlabel(handles.axes4, 'Time (s)');
        ylabel(handles.axes4, 'Deflection (mm)');
        legend(handles.axes4, 'Experiment', 'Simulation (reset)', 'Location', 'northeast');
        
        box(handles.axes4, 'on');
        xlim(handles.axes4, [t_min, t_max]); 
        set(handles.axes4, 'FontSize', 10);
        set(handles.axes4, 'FontName', 'Arial');
        
               
        
end