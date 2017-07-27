
function faml_plot_figure_reset(legend_m, data, data_calibrated, t_min, t_max)
        
        dim = size(data_calibrated);
        num_sim = dim(2);
    
    
        h1 = figure(20);
        h2 = figure(21);
        h3 = figure(22);
                
        figure(h1);
        cla(h1);
        
        plot(data(:,1), data(:,2), 'b');        
        hold on;
        
        for i = 1:num_sim
            
            plot(data_calibrated{i}(:,1), data_calibrated{i}(:,2), 'color', rand(1,3));
            hold on;
        
        end
         
        title('Force vs Time (Exp vs Sim (reset))');
        xlabel('Time (sec)');
        ylabel('Force (kN)');
        legend(legend_m, 'Location', 'northeast');
         
        box on; 
        %xlim([t_min, t_max]); 
        set(gca, 'FontSize', 15);
        set(gca, 'FontName', 'Arial');
                
                
        figure(h2);
        cla(h2);
        
        plot(data(:,1), data(:,3), 'b');        
        hold on;
        
        for i = 1:num_sim
            
            plot(data_calibrated{i}(:,1), data_calibrated{i}(:,3), 'color', rand(1,3));
            hold on;
            
        end
         
        title('Displacement vs Time (Exp vs Sim (reset))');
        xlabel('Time (sec)');
        ylabel('Displacement (mm)');
        legend(legend_m, 'Location', 'northeast');
         
        box on; 
        %xlim([t_min, t_max]); 
        set(gca, 'FontSize', 15);
        set(gca, 'FontName', 'Arial');
                
          
        figure(h3);
        cla(h2);
        
        plot(data(:,2), data(:,3), 'b');        
        hold on;
        
        for i = 1:num_sim
            
            plot(data_calibrated{i}(:,2), data_calibrated{i}(:,3), 'color', rand(1,3));
            
        end
         
        title('Displacement vs Force (Exp vs Sim (reset))');
        xlabel('Force (kN)');
        ylabel('Displacement (mm)');
        legend(legend_m, 'Location', 'northeast');
         
        box on; 
        %xlim('auto'); 
        set(gca, 'FontSize', 15);
        set(gca, 'FontName', 'Arial');
                       
end 
        
