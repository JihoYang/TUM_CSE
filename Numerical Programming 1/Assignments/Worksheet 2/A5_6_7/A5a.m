function result=A5a(c,v)
% naïve: create circulant matrix and multiply with v
m=numel(c);
mm1 = m-1;
C = c(mod( bsxfun(@minus,(0:mm1)',0:mm1), m )+1);
result = C*v;
end


function result2 = DFT(c,v)

result2 = ifft(fft(c)*fft(v));



end