
%FAML_PLOT_GUI generates plots for imported data at the pre-defined
%axes in the GUI environment

function faml_plot_GUI(time, data, t_min, t_max, data_type, handles)
    
    
    if data_type == 'Experiment';
    
        plot(handles.axes1, time(:), data(:,4));
        title(handles.axes1, '4th Sensor Disp (-z) vs Time');
        xlabel(handles.axes1, 'Time (s)');
        ylabel(handles.axes1, '4th Sensor Displacement (mm)');
        legend(handles.axes1, 'Disp (4th) vs Time');
        box(handles.axes1, 'on');
        xlim(handles.axes1, [t_min, t_max]);
        set(handles.axes1, 'FontSize', 10);
        set(handles.axes1, 'FontName', 'Arial');
  
       
        plot(handles.axes2, time(:), data(:,7));
        title(handles.axes2, '5th Sensor Accel (+z) vs Time');
        xlabel(handles.axes2, 'Time (s)');
        ylabel(handles.axes2, '5th Sensor Accleration (mm/s^2)');
        legend(handles.axes2, 'Accel (5th) vs Time');
        box(handles.axes2, 'on');
        xlim(handles.axes2, [t_min, t_max]);
        set(handles.axes2, 'FontSize', 10);
        set(handles.axes2, 'FontName', 'Arial');
        
        
        plot(handles.axes3, time(:), data(:,10));
        title(handles.axes3, '6th Sensor Accel (+z) vs Time');
        xlabel(handles.axes3, 'Time (s)');
        ylabel(handles.axes3, '6th Sensor Accleration (mm/s^2)');
        legend(handles.axes3, 'Accel (6th) vs Time');
        box(handles.axes3, 'on');
        xlim(handles.axes3, [t_min, t_max]);
        set(handles.axes3, 'FontSize', 10);
        set(handles.axes3, 'FontName', 'Arial');
        
    else
        
        plot(handles.axes4, time, data(:,2), 'r');
        title(handles.axes4, 'Deflection comparison (Exp vs Sim)');
        legend(handles.axes4, 'Experiment', 'Simulation');
        xlim(handles.axes4, [t_min, t_max]);
        
        
    end
        
end






