
%plot_figure function creates external plot figures for imported data

function faml_plot_figure(time, data, defl_exp, t_min, t_max, data_type)
        
    if data_type == 'Experiment'
    
        figure(1);
        plot(time, data(:,4));
        title('4th Sensor Disp (-z) vs Time');
        xlabel('Time (s)');
        ylabel('4th Sensor Displacement (mm)');
        legend('Disp (4th) vs Time', 'Location', 'northeast');
        set(gca, 'FontSize', 15);
        set(gca, 'FontName', 'Arial');
        xlim([t_min, t_max]);
        box on;
        
        
        figure(2);       
        plot(time, data(:,7));
        title('5th Sensor Accel (+z) vs Time');
        xlabel('Time (s)');
        ylabel('5th Sensor Accleration (mm/s^2)');
        legend('Accel (5th) vs Time', 'Location', 'northeast');
        set(gca, 'FontSize', 15);
        set(gca, 'FontName', 'Arial');
        xlim([t_min, t_max]);
        box on;
        
        
        figure(3);       
        plot(time, data(:,10));
        title('6th Sensor Accel (+z) vs Time');
        xlabel('Time (s)');
        ylabel('6th Sensor Accleration (mm/s^2)');
        legend('Accel (6th) vs Time', 'Location', 'northeast');
        set(gca, 'FontSize', 15);
        set(gca, 'FontName', 'Arial');
        xlim([t_min, t_max]);
        box on;
        
        
    else
        
        figure(4);
        h = figure(4);
        
        if isempty(get(h, 'Children'))
            
            plot(defl_exp(:,1), defl_exp(:,2), 'b');
            hold on;
            plot(time, data(:,2), 'r');
            
            title('Deflection comparison (Exp vs Sim)');
            xlabel('Time (s)');
            ylabel('Deflection (mm)');
            legend('Experiment', 'Simulation');
            xlim([t_min, t_max]);
            
        else
           
            plot(time, data(:,2), 'r');
            title('Deflection comparison (Exp vs Sim)');
            xlabel('Time (s)');
            ylabel('Deflection (mm)');
            legend('Experiment', 'Simulation');
            xlim([t_min, t_max]);
            
        end
        
   
    end
    
        

end









