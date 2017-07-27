function A21_showData (filename,E,C)
% shows measurements (a) and K(T) of (d) if E and C are given
% function A21_showData (filename,E,C)
% Input
%   filename        filename with measurements
%   E,C             (optional) parameters for graph K(T)


M=load(filename);                    % load measurements
if nargin==3                         % E,C are given => T range for graph
    T=linspace(min(M(:,1))-1,max(M(:,1))+1,500);
end

subplot(2,1,1);
semilogy(1./M(:,1),M(:,2),'k.');     % plot measurements as dots
set(gca,'YScale','log');
if nargin==3                         % if C and E are given 
    hold on;semilogy(1./T,C*exp(-E./T));hold off;
    title(sprintf('%s, E=%g, C=%g',filename,E,C));
end
grid on; xlabel('1/T');ylabel('K');

subplot(2,1,2);
plot(M(:,1),M(:,2),'k.');
if nargin==3                         % if C and E are given 
    hold on;plot(T,C*exp(-E./T));hold off;
end
xlabel('T');ylabel('K');

end
