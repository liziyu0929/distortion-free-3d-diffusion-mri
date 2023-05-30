function [res,RESVEC] = pcgSPIRiTL1_phaseCorrection(kdata, phs, sampling_order, GOP, XOP, nIter,nIterSplit, sz, nch, z_start, lambda,alpha,beta,x0)

imSize = size(x0);

% make dyadic size if Wavelets are used. 
if strcmp(class(XOP),'Wavelet') == 1

    if length(imSize)>2
        imSize_dyd = [max(2.^ceil(log2(imSize(1:2)))), max(2.^ceil(log2(imSize(1:2)))),imSize(3)];
    else
        imSize_dyd = [max(2.^ceil(log2(imSize(1:2)))), max(2.^ceil(log2(imSize(1:2))))];
    end
    
else
    imSize_dyd = imSize;
end

y = epi2D_phs_fw_adj(kdata,phs,sampling_order,z_start);
y = fft2c(y);

N = prod(sz)*nch;
z = zeros([sz nch]);
tmpres = [x0(:); z(:)];

for n = 1:nIterSplit
    b = [y(:); sqrt(beta)*tmpres(1:N)];
    
    [tmpres,FLAG,RELRES,ITER,RESVEC] = pcg(@(x)afun(x,phs,sampling_order,GOP,sz,nch,z_start,beta,alpha),...
        b,1e-6,round(nIter/nIterSplit),[],[],tmpres(:));
    
    res = reshape(tmpres(1:N),[sz nch]);
    resL1 = ifft2c(res);
    tmp = zpad(resL1,imSize_dyd);
    tmp = XOP*tmp;
    tmp = SoftThresh(tmp,lambda/sqrt(beta));
    resL1 = XOP'*tmp;
    resL1 = reshape(resL1,imSize_dyd);
    resL1 = fft2c(crop(resL1,imSize));
    
    tmpres = [resL1(:); tmpres(N+1:end)];
    
end


    function y = afun(x,phs,sampling_order,GOP,sz,nch,z_start,beta,alpha)
        
        l = prod(sz)*nch;
        x_k_in = reshape(x(1:l), [sz nch]);
        x_im = ifft2c(x_k_in);
        x_k = epi2D_phs_fw(x_im, phs, sampling_order,z_start); % DP
        x_im = epi2D_phs_fw_adj(x_k,phs,sampling_order,z_start); % (DP)^H
        x_k = fft2c(x_im);
        x_calib = GOP'*(GOP*x_k_in);
        y = [x_k(:)+alpha*x_calib(:); sqrt(beta)*x(1:l)];
    end
end