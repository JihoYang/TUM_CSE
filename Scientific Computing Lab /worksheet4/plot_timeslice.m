
function plot_timeslice(mesh,n,t_length,T, row, column,fig, sub_title)
% Function which plots the whole timeslice figure

% m : mesh [X,Y] for all N. Double cell array.
% n : Number of points in grid;
%  t: timeslice
%  T: solutions for this t
% [row,column]: size of subplot;
% fig : Figure number;
% sub_title : Title of figure.

%PlotTitle = [result.method, ' at ', sub_title];
%hold on;

figure(fig);
hold on;
s = 1;
N = 0;
for i = 1:t_length

    % Get the correct mesh size
    if N  ~= n(i)
        N = n(i);
        X = mesh{s}{1};
        Y = mesh{s}{2};
        s = s+1;
    end

    subplot(row,column, i);
    surf(X,Y,T{i});

    if(i ==1)
        xlabel('x');
        ylabel('y');
        zlabel('T');
    end
    %axis([0,1,0,1,0,1]);
end

set(gcf,'position',get(0,'screensize'))


subtitle(sub_title);

%Save plots

%Change the incompatible characters into underbar
filename = char(strrep(sub_title(:,:), ',', ''));
%filename = char(strrep(filename(:,:), '.', '_'));
filename = char(strrep(filename(:,:), ':', ''));

saveas(gcf, [pwd '/figures/', filename], 'jpeg');
