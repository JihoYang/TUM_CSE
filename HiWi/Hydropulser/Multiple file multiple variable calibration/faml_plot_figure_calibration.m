
function faml_plot_figure_calibration(data, data_calibrated, ptitle, xoutput, youtput, t_min, t_max, handles)
        
    if get(handles.calibrate_force, 'Value')
       
       k = 1;
     
    end
    
    if get(handles.calibrate_disp, 'Value')
       
       k = 2;
     
    end
    
    if get(handles.calibrate_disp_force, 'Value')
       
       k = 3;
     
    end
    
    h = figure(k+19);
    cla(h);
        
    figure(k+19);
        
    plot(data(:,1), data(:, k+1), 'b');
    hold on;
    plot(data_calibrated(:,1), data_calibrated(:,2), 'r');
            
    title(ptitle(k));
    xlabel(xoutput(k));
    ylabel(youtput(k));
    legend('Experiment', 'Simulation (Calibrated)', 'Location', 'northeast');
        
    box on; 
    xlim([t_min, t_max]);
    
    set(gca, 'FontSize', 15);
    set(gca, 'FontName', 'Arial');
    
    if k == 3
            
        h = figure(k+19);
        cla(h);
        
        figure(k+19);
        
        plot(data(:,2), data(:, 3), 'b');
        hold on;
        plot(data_calibrated(:,1), data_calibrated(:,2), 'r');
            
        title(ptitle(k));
        xlabel(xoutput(k));
        ylabel(youtput(k));
        legend('Experiment', 'Simulation (Calibrated)', 'Location', 'northeast');
        
        box on; 
        xlim('auto');
    
        set(gca, 'FontSize', 15);
        set(gca, 'FontName', 'Arial');
           
    end
        
end
