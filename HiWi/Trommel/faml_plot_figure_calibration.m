
function faml_plot_figure_calibration(time, data_exp, data_calibrated, t_min, t_max, type, legend_m, handles)
    
    if type == 'Experiment'
        
        cla(figure(4));
        figure(4);
        plot(time, data_calibrated(:,2), 'b');
        title('Deflection vs Time');
        xlabel('Time (s)');
        ylabel('Deflection (mm)');
        legend('Deflect vs t', 'Location', 'northeast');
        xlim([t_min, t_max]);
        set(gca, 'FontSize', 15);
        set(gca, 'FontName', 'Arial');
        box on;
        hold on;
        
    else
        
        cla(figure(4));
        figure(4);
        plot(data_exp(:,1), data_exp(:,2), 'b');
        hold on;
        plot(data_calibrated(:,1), data_calibrated(:,2), 'r');
        title('Deflection Comparison');
        xlabel('Time (s)');
        ylabel('Deflection (mm)');
        legend(legend_m);
        xlim([t_min, t_max]);
        set(gca, 'FontSize', 15);
        set(gca, 'Fontname', 'Arial');
        box on;
        
end