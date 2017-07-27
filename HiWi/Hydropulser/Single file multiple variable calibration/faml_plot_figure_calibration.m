
function faml_plot_figure_calibration(legend_m, data, data_calibrated, ptitle, xoutput, youtput, t_min, t_max, handles)
    
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
        
    h = figure(k+19);
    cla(h);
        
    figure(k+19);
        
    plot(data(:,1), data(:, k+1), 'b');
    hold on;
    
    for i = 1:num_sim
        
        plot(data_calibrated{i}(:,1), data_calibrated{i}(:,2), 'color', rand(1,3));
        hold on;
        
    end
            
    title(ptitle(k))
    xlabel(xoutput(k));
    ylabel(youtput(k));
    legend(legend_m);
        
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
        
        for i = 1 : num_sim
            
            plot(data_calibrated{i}(:,1), data_calibrated{i}(:,2), 'color', rand(1,3));
            hold on;
        
        end
            
        title(ptitle(k));
        xlabel(xoutput(k));
        ylabel(youtput(k));
        legend(legend_m, 'Location', 'northeast');
        
        box on;
        xlim('auto');
    
        set(gca, 'FontSize', 15);
        set(gca, 'FontName', 'Arial');
        hold on;
           
    end
    
    
        
end