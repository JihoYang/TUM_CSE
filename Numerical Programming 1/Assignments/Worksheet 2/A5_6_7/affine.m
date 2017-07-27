function y = affine(a,b,n)
    
    
    
    for i = 1:n
        
        x(i) = cos(pi*i/n);
        y(i) = a+(b-a)/2+(b-a)/2*x(i);
        
    end
    
    
        
        


end