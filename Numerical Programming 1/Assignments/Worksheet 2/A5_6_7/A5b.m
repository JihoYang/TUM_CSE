function z=A5b(c,v)
% calculate C*v using fft and ifft.
tic;
z = ifft( fft(c).*fft(v) );
toc;
disp(toc);
end