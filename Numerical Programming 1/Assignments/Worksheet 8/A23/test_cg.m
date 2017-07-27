function [x_iter,omega]=test_cg(A,b,max_iter)
normA = norm(A,inf); 
x_iter=zeros(numel(b),max_iter); omega=zeros(1,max_iter);

cg(@(x) A*x,b,zeros(size(b)),@check_tol);

    function res=check_tol(x,res,k)
        if k>0
            x_iter(:,k)=x;
            omega(k)=norm(res,inf)/(normA*norm(x,inf));
        end
        res=k>=max_iter;
    end

end