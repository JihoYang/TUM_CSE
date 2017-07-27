
function faml_plot_figure_reset(data, data_calibrated, t_min, t_max)
        
        h1 = figure(20);
        h2 = figure(21);
        h3 = figure(22);
                
        figure(h1);
        cla(h1);
        
        plot(data(:,1), data(:,2), 'b');        
        hold on;
        plot(data_calibrated(:,1), data_calibrated(:,2), 'r');
         
        title('Force vs Time (Exp vs Sim (reset))');
        xlabel('Time (sec)');
        ylabel('Force (kN)');
        legend('Experiment', 'Simulation (original)', 'Location', 'northeast');
         
        box on; 
        xlim([t_min, t_max]); 
        set(gca, 'FontSize', 15);
        set(gca, 'FontName', 'Arial');
                
                
        figure(h2);
        cla(h2);
        
        plot(data(:,1), data(:,3), 'b');        
        hold on;
        plot(data_calibrated(:,1), data_calibrated(:,3), 'r');
         
        title('Displacement vs Time (Exp vs Sim (reset))');
        xlabel('Time (sec)');
        ylabel('Displacement (mm)');
        legend('Experiment', 'Simulation (original)', 'Location', 'northeast');
         
        box on; 
        xlim([t_min, t_max]); 
        set(gca, 'FontSize', 15);
        set(gca, 'FontName', 'Arial');
                
          
        figure(h3);
        cla(h2);
        
        plot(data(:,2), data(:,3), 'b');        
        hold on;
        plot(data_calibrated(:,2), data_calibrated(:,3), 'r');
         
        title('Displacement vs Force (Exp vs Sim (reset))');
        xlabel('Force (kN)');
        ylabel('Displacement (mm)');
        legend('Experiment', 'Simulation (original)', 'Location', 'northeast');
         
        box on; 
        xlim('auto'); 
        set(gca, 'FontSize', 15);
        set(gca, 'FontName', 'Arial');
                       
end 
        
