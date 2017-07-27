function A4
a_unstable = @(x) log(1+x);

%% Plot 1
figure(1);clf;
x=(1:500)/50*eps;
plot(x,a_unstable(x),x,a_stable(x));
legend({'unstable','stable'},'Location','NorthWest');

%% Plot 2
figure(2);clf;
x=(1:5000)/500*eps;
err = @(alg,true) max(abs((alg-true)./(true)),eps);
semilogy(x,err(a_unstable(x),log1p(x)),x,err(a_stable(x),log1p(x)));
grid on; title('relative error');
legend({'unstable','stable'});

end

function erg=a_stable(x)
w=1+x;
erg=(log(w)./(w-1)).*x;
erg(w==1)=x(w==1);
end
