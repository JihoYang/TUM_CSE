
function faml_plot_GUI_calibration(time, data_exp, data_calibrated, t_min, t_max, type, legend_m, handles)
    
    if type == 'Experiment'
        
        cla(handles.axes4);
        plot(handles.axes4, time, data_calibrated(:,2), 'b'); 
        title(handles.axes4, 'Deflection vs Time');
        xlabel(handles.axes4, 'Time (s)');
        ylabel(handles.axes4, 'Deflection (mm)');
        legend(handles.axes4, 'Deflect vs t', 'Location', 'northeast');   
        box(handles.axes4, 'on');
        xlim(handles.axes4, [t_min, t_max]); 
        set(handles.axes4, 'FontSize', 10);
        set(handles.axes4, 'FontName', 'Arial');
        hold(handles.axes4, 'on');
        
    else
        
        cla(handles.axes4);
        plot(handles.axes4, data_exp(:,1), data_exp(:,2), 'b');
        hold(handles.axes4, 'on');
        plot(handles.axes4, data_calibrated(:,1), data_calibrated(:,2), 'r');
        title(handles.axes4, 'Deflection Comparison');
        xlabel(handles.axes4, 'Time (s)');
        ylabel(handles.axes4, 'Deflection (mm)');
        legend(handles.axes4, legend_m);
        box(handles.axes4, 'on');
        xlim(handles.axes4, [t_min, t_max]);
        set(handles.axes4, 'FontSize', 10);
        set(handles.axes4, 'FontName', 'Arial');
        hold(handles.axes4, 'on');
        
    end
            
        
end