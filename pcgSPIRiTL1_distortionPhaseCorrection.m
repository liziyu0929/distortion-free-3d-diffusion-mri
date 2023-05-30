function [res,RESVEC] = pcgSPIRiTL1_distortionPhaseCorrection(kdata_ap, kdata_pa,fm,t_es, phs_ap,phs_pa,samp_ap, samp_pa, GOP, XOP, nIter,nIterSplit, sz, nch,lambda,alpha,beta,x0)

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

y = epi2D_dis_phs_fw2_adj(kdata_ap,kdata_pa,t_es,fm,phs_ap, phs_pa, samp_ap,samp_pa);
y = fft2c(y);

N = prod(sz)*nch;
z = zeros([sz nch]);
tmpres = [x0(:); z(:)];

for n = 1:nIterSplit
    b = [y(:); sqrt(beta)*tmpres(1:N)];
    
    [tmpres,FLAG,RELRES,ITER,RESVEC] = pcg(@(x)afun(x,fm,t_es, phs_ap,phs_pa,samp_ap, samp_pa,GOP,sz,nch,beta,alpha),...
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


    function y = afun(x,fm,t_es, phs_ap,phs_pa,samp_ap, samp_pa,GOP,sz,nch,beta,alpha)
        
        l = prod(sz)*nch;
        x_k = reshape(x(1:l), [sz nch]);
        x_im = ifft2c(x_k);
        
        kdata_ap_tmp = epi2D_dis_phs_fw2(x_im,t_es,fm,phs_ap,samp_ap,1);
        kdata_pa_tmp = epi2D_dis_phs_fw2(x_im,t_es,fm,phs_pa,samp_pa,0);
        
        yy = epi2D_dis_phs_fw2_adj(kdata_ap_tmp,kdata_pa_tmp,t_es,fm,phs_ap, phs_pa, samp_ap,samp_pa);
        
        yy = fft2c(yy);
        
        x_calib = GOP'*(GOP*x_k);
        
        y = [yy(:)+alpha*x_calib(:); sqrt(beta)*x(1:l)];
        
    end
end