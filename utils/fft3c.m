function b=fft3c(a)


b=fft2c(a);
b=fft1c(b,3);

end