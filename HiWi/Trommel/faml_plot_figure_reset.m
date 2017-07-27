
function faml_plot_figure_reset(data, data_calibrated, t_min, t_max)
        
        h1 = figure(4);
       
        figure(h1);
        cla(h1);
        
        plot(data(:,1), data(:,2), 'b');        
        hold on;
        plot(data_calibrated(:,1), data_calibrated(:,2), 'r');
         
        title('Deflection Comparison (Exp vs Sim (reset))');
        xlabel('Time (s)');
        ylabel('Deflection (mm)');
        legend('Experiment', 'Simulation (original)', 'Location', 'northeast');
         
        box on;
        xlim([t_min, t_max]); 
        set(gca, 'FontSize', 15);
        set(gca, 'FontName', 'Arial');
                
                       
end 
        
