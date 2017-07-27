


% Let's take a closer look at accuracy vs storage, runtime

inputs = table(granularities', mesh, forcing, boundary);
inputs.Properties.VariableNames{1} = 'size';

rows = height(performance);
[error, reduction] = deal(zeros(rows,1));
for row=1:rows
    
    index = find(inputs.size==performance.size(row));
    
    error(row) = get_error(performance.solution{row},exact_soln{index});
    
    if (index ~= 1)
        reduction(row) = error(row-1)/error(row);
    end
end
outputs = table(performance.method, performance.size, performance.runtime, performance.storage, error, reduction);
outputs.Properties.VariableNames = {'method', 'size', 'runtime', 'storage', 'error', 'reduction'};
outputs

figure
methods = unique(outputs.method);
for m=1:length(methods)
    grouping = outputs(strcmp(outputs.method,methods(m)),:);
    loglog(grouping.error, grouping.runtime, '*-')
    hold on;
end
legend(methods);
grid on;
xlabel 'Error'
ylabel 'Runtime'


figure
methods = unique(outputs.method);
for m=1:length(methods)
    grouping = outputs(strcmp(outputs.method,methods(m)),:);
    loglog(grouping.error, grouping.storage, '*-')
    hold on;
end
legend(methods);
grid on;
xlabel 'Error'
ylabel 'Storage'

