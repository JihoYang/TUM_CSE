load('Adams_result_single_body.tab');
time_ref = Adams_result_single_body(:,1);
time_ref = time_ref(1:100);
disp_ref = Adams_result_single_body(:,2);
disp_ref = disp_ref(1:100);
load('out_Euler.dat');
time_Euler = out_Euler(:,1);
disp_Euler = out_Euler(:,2);
disp_Euler(1000) = disp_Euler(999);
disp_Euler = disp_Euler(1:10:1000);
load('out_RK.dat');
time_RK = out_RK(:,1);
disp_RK = out_RK(:,2);
disp_RK(1000) = disp_RK(999);
disp_RK = disp_RK(1:10:1000);
load('out_BDF.dat');
time_BDF = out_BDF(:,1);
disp_BDF = out_BDF(:,2);
disp_BDF(1000) = disp_BDF(999);
disp_BDF = disp_BDF(1:10:1000);

plot(time_ref,disp_ref,'r');
hold on;
plot(time_ref,disp_Euler,'b');
plot(time_ref,disp_RK,'k');
plot(time_ref,disp_BDF,'g');

h = legend('Reference from ADAMS','Explicit Euler','Runge Kutta','Backward Differentiation Formula/GEAR');
set(gca,'FontSize',14)
xlabel('Time', 'FontSize',14);
ylabel('Displacement in mm', 'FontSize', 14);
resize_legend(h, 1.1);
fn = 'Comparission_of_Solvers';
saveas(gcf, fn, 'fig');
print( gcf, '-dbmp16m', fn );

