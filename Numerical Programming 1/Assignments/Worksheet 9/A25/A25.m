function A25

x0=1/sqrt(3);

[z,fail,steps]=NewtonRelax(@F_func,...
    @(z) adiffget(F_func(adiff(z)),'derivatives'),...
    [x0;x0],1e-15);

if fail
    error('NewtonRelax did not converge');
end

display(z);
display(steps);



end

function J=DF_func(z)
J=adiffget(F_func( adiff(z) ),'derivatives');
end