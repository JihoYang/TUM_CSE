function A5c

sizerange = 2.^(8:13);
comptime = zeros(2,numel(sizerange));

for k=1:numel(sizerange)
    s = sizerange(k);
    c = rand(s,1); v=rand(s,1);
    tic; A5a(c,v);  comptime(1,k)=toc;
    tic; A5b(c,v);  comptime(2,k)=toc;
end

sprintf(' m =%5i   time_a = %1.5f  time_b = %1.5f\n',[sizerange;comptime])

end