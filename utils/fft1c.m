function Y = fft1c(X,dim)

if nargin < 2
    dim=1;
end

l = size(X,dim);

Y = ifftshift(fft(fftshift(X,dim),[],dim),dim)/sqrt(l);

end