function b=ifft3c(a)

b=ifft2c(a);
b=ifft1c(b,3);


end