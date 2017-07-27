
%plot_figure function creates external plot figures for imported data

function faml_plot_figure(legend_m_sim, data, t_min, t_max, variables, unit, data_type)
        
    if data_type == 'Experiment'
        
        for k = 1:2
            
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
            
        figure(4);
        plot(data(:,1), data(:,2));
        title(strcat(variables(1), ' variation over time'));
        xlabel('Time (sec)');
        ylabel(strcat(variables(1), unit(1)));
        legend(legend_m_sim{:,1}, 'Location', 'northeast');
        set(gca, 'FontSize', 15);
        set(gca, 'FontName', 'Arial');
        xlim([t_min, t_max]);
        box on; 
        hold on;
        
        
        figure(5);
        plot(data(:,1), data(:,3));
        title(strcat(variables(2), ' variation over time'));
        xlabel('Time (sec)');
        ylabel(strcat(variables(2), unit(2)));
        legend(legend_m_sim{:,2}, 'Location', 'northeast');
        set(gca, 'FontSize', 15);
        set(gca, 'FontName', 'Arial');
        xlim([t_min, t_max]);
        box on; 
        hold on;
            
    
    figure(6);
        
    plot(data(:,2), data(:,3));
        
    title('Displacement variation over Force');
    xlabel('Force (kN)');
    ylabel('Displacement (mm)');
    legend(legend_m_sim{:,3}, 'Location', 'northeast');
    set(gca, 'FontSize', 15);
    set(gca, 'FontName', 'Arial');
        
    box on; 
    hold on;
    
    figure(7);
    plot(data(:,1), data(:,4));
    legend(legend_m_sim{:,4});
    box on; 
    hold on;
    
    figure(8);
    plot(data(:,1), data(:,5));
    legend(legend_m_sim{:,5});
    box on;
    hold on;
        
    end
        

end









