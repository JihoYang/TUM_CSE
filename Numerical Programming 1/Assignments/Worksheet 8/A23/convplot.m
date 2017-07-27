function convplot(A,b,max_iter)
% Plot convergence plots for (relaxed) Jacobi/SOR and CG method.

names = {'Jacobi λ=0.5', 'Jacobi λ=1.0', 'Jacobi λ=1.5',...
         'SOR λ=0.5', 'SOR λ=1.0', 'SOR λ=1.5', 'CG'};
funcs = { @() classicIteration('Jacobi',A,b,0.5,max_iter), ...
          @() classicIteration('Jacobi',A,b,1.0,max_iter), ...
          @() classicIteration('Jacobi',A,b,1.5,max_iter), ...
          @() classicIteration('SOR'   ,A,b,0.5,max_iter), ...
          @() classicIteration('SOR'   ,A,b,1.0,max_iter), ...
          @() classicIteration('SOR'   ,A,b,1.5,max_iter), ...
          @() test_cg(A,b,max_iter) }; % set of algorithms
anz=numel(names);
x_direct = A\b; n = size(A,1);         % direct solution; for forward error
      
err_backward=NaN*ones(anz,max_iter);
err_forward=err_backward;
for k=1:anz
    [x_iter,omega]=funcs{k}();
    err_backward(k,:)=omega;
    err_forward(k,:)=max(abs(x_iter-repmat(x_direct,1,max_iter)),[],1);
end
err_backward(err_backward<eps)=eps;
err_forward(err_forward<eps)=eps;

subplot(2,1,1);
semilogy(1:max_iter,err_backward);grid on;legend(names);
title(['backward error; n=',num2str(n),...
    ' density=',num2str(nnz(A)/(n*n),2),...
    ' cond≈',num2str(condest(A),2)]);
          
subplot(2,1,2);
semilogy(1:max_iter,err_forward);grid on;
axis([1,max_iter,1e-16,1e10]);
title('forward error');

end