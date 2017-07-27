    [filename_exp pathname_exp] = uigetfile({'.asc'},'File Selector');
    fullpathname_exp = strcat(pathname_exp, filename_exp);

    %Build full file name and read/import it to MATLAB GUI workspace
    loaddata_exp = fullfile(pathname_exp, filename_exp);
    data_cellarray = textread(loaddata_exp,'%s','delimiter','\n');
    
    %Remove header (assuming the first row of data starts from the 39th)
    data_cellarray(1:38) = [];
    
    %Change commas to dots and create numeric matrix
    data_numeric_exp = double(str2num(char(strrep(data_cellarray,',','.')')));
    
    
    
    %Define/Compute parameters for data plotting
    freq = 100;
    dimension_exp = size(data_numeric_exp);
    t_exp = zeros(dimension_exp(1),1);
    t_exp(2:end) = 1/freq*(1:dimension_exp(1)-1);      
    
    %offset_value = 7.771;
    %correction_factor = 1.25;
    
    offset_value = 0;
    correction_factor= 1;
    
    data_exp = data_numeric_exp;
    data_exp(:,4) = (data_exp(:,4)+offset_value).*correction_factor;
    
    defl_exp(:,1) = t_exp;
    defl_exp(:,2) = data_exp(:,4);
    
    [pks loc] = findpeaks(defl_exp(:,2));
    number_peak = size(pks);
    avg_peak = mean(pks);

    while number_peak(1) > 10
        
        pks = findpeaks(pks);
        number_peak = size(pks);
        avg_peak = mean(pks);
        
    end
    
pks_max = max(pks);
pks_min = min(pks);
pks_diff = abs(pks_max - pks_min);

remove = find(pks < pks_max - 0.9*pks_diff);
pks(remove) = [];
number_peak = size(pks);

%Remove initial oscillation before first peak
first_peak = find(defl_exp(:,2) == pks(1));
defl_exp((1:first_peak-1),:) = [];

%Compute excitation period
first_peak = find(defl_exp(:,2) == pks(1));
second_peak = find(defl_exp(:,2) == pks(2));
period = t_exp(second_peak) - t_exp(first_peak);
number_excitation = floor((defl_exp(end,1)-defl_exp(1,1))/period);

cycle = find(abs(defl_exp(:,1) - (defl_exp(1,1) + period)) < 1000*eps);

    for i = 1 : cycle - 1
        
        for j = 1 : number_excitation
            
            excit(i,j) = defl_exp((j-1)*(cycle-1)+i,2);
            
        end
        
    end
    
%Calculate the average excitation cycle
excit_tp = excit';
avg_defl_y = mean(excit_tp)';
avg_defl_x = defl_exp((1:cycle-1),1);
avg_defl(:,1) = avg_defl_x;
avg_defl(:,2) = avg_defl_y;

%Remove last bit of uncessary data
avg_defl((end-20:end),:) = [];

%Shift the averaged deflection's time so it starts from t = 0s
min_time_avg_defl = min(avg_defl(:,1));
avg_defl(:,1) = avg_defl(:,1) - min_time_avg_defl;

%damping coefficient calculation
peak_max = max(avg_defl(:,2));
locs_max = find(peak_max - avg_defl(:,2) < 1000*eps);
[peak_avg locs_avg] = findpeaks(avg_defl(:,2));
dim_peak_avg = size(peak_avg);
peak = zeros(dim_peak_avg(1)+1,1);
locs = zeros(dim_peak_avg(1)+1,1);
peak(1) = peak_max;
peak(2:end) = peak_avg;
locs(1) = locs_max;
locs(2:end) = locs_avg;
P = polyfit(avg_defl(locs,1), log(peak),1);
delta = -P(1);
%mass = force/9.81;
%damping_coefficient = 2*delta*mass;

%Damping ratio calculation - OPTIONAL
%dim_avg_defl = size(avg_defl);
%dim_peak = size(peak);
%cycle_damping = dim_peak(1);
%zeta = (log(peak(1)/peak(end)))/((cycle_damping-1)*2*pi);

plot(avg_defl(locs,1), peak, 'ro');
 hold on;
 plot(avg_defl(:,1), avg_defl(:,2));
 plot(avg_defl(:,1), exp(P(2)+P(1).*avg_defl(:,1)), 'color', [0 0.5 0]);
 title('Exponentially decaying vibration');                                                                                                                                       
 xlabel('Time (s)'); 
 ylabel('Deflection (mm)');
 box on;
 legend('Peak', 'Averaged deflection', 'Exponential curve', ...
         'Location', 'northeast');
 set(gca, 'FontName', 'Arial');
 set(gca, 'FontSize', 15);
