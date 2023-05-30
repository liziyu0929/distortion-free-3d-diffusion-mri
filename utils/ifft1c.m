function Y = ifft1c(X,dim)

if nargin<2
    dim=1;
end
l = size(X,dim);

Y = fftshift(ifft(ifftshift(X,dim),[],dim),dim)*sqrt(l);



end