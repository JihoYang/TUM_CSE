
%plot_figure function creates external plot figures for imported data

function faml_plot_figure(data, dimension, t_min, t_max, variables, unit, data_type)
        
    if data_type == 'Experiment'
        
        for k = 1:dimension(2)-2
            
            figure(k);
        
            plot(data(:,1), data(:,k+1));
        
            title(strcat(variables(k), ' variation over time'));
            xlabel('Time (sec)');
            ylabel(strcat(variables(k), unit(k)));
            legend(strcat(variables(k), ' vs Time'), 'Location', 'northeast');
            set(gca, 'FontSize', 15);
            set(gca, 'FontName', 'Arial');
            xlim([t_min, t_max]);
        
            box on; 
        
        end
    
    figure(3);
    
    plot(data(:,2), data(:,3));
            
    title('Displacement variation over Force');
    xlabel('Force (kN)');
    ylabel('Displacement (mm)');
    legend('Displacement vs Force', 'Location', 'northeast');
    set(gca, 'FontSize', 15);
    set(gca, 'FontName', 'Arial');
    
    box on;
        
    else
        
        for k = 1:2
            
            figure(k+3);
        
            plot(data(:,1), data(:,k+1));
        
            title(strcat(variables(k), ' variation over time'));
            xlabel('Time (sec)');
            ylabel(strcat(variables(k), unit(k)));
            legend(strcat(variables(k), ' vs Time'), 'Location', 'northeast');
            set(gca, 'FontSize', 15);
            set(gca, 'FontName', 'Arial');
            xlim([t_min, t_max]);
        
            box on; 
        
        end
    
    figure(6);
        
    plot(data(:,2), data(:,3));
        
    title('Displacement variation over Force');
    xlabel('Force (kN)');
    ylabel('Displacement');
    legend('Displacement vs Force', 'Location', 'northeast');
    set(gca, 'FontSize', 15);
    set(gca, 'FontName', 'Arial');
        
    box on; 
        
    end
        

end









