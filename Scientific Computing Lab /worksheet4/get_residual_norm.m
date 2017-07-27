
function R = get_residual_norm(b,T,wx,wy,wh)

    dim = size(T);
    ny = dim(1) -2;
    nx = dim(2) -2;
    R_mat = zeros(dim);

    for r = 2:dim(1) -1
        for c = 2:dim(2)-1
            R_mat(r,c)= -b(r,c) + wh*T(r,c)...
            +wx*(T(r,c-1)+T(r,c+1))...
            +wy*(T(r-1,c)+T(r+1,c));
        end
    end

    R = sqrt((1/(nx*ny))*sum(sum(R_mat.^2)));
        
end

