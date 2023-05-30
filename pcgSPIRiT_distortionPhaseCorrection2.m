function [res,RESVEC] = pcgSPIRiT_distortionPhaseCorrection2(kdata_ap, kdata_pa,fm,t_es, phs_ap,phs_pa,samp_ap, samp_pa, GOP, nIter, sz, nch,lambda,alpha)

y = epi2D_dis_phs_fw2_adj(kdata_ap,kdata_pa,t_es,fm,phs_ap, phs_pa, samp_ap,samp_pa);
y = fft2c(y);

N = prod(sz)*nch;
z = zeros([sz nch]); % Tikhonov regularization term
b = [y(:); z(:)];

[tmpres,FLAG,RELRES,ITER,RESVEC] = pcg(@(x)afun(x,fm,t_es, phs_ap,phs_pa,samp_ap, samp_pa,GOP,sz,nch,lambda,alpha),b,1e-6,nIter);


res = reshape(tmpres(1:N),[sz nch]);

    function y = afun(x,fm,t_es, phs_ap,phs_pa,samp_ap, samp_pa,GOP,sz,nch,lambda,alpha)
        
        l = prod(sz)*nch;
        x_k = reshape(x(1:l), [sz nch]);
        x_im = ifft2c(x_k);
        
        kdata_ap_tmp = epi2D_dis_phs_fw2(x_im,t_es,fm,phs_ap,samp_ap,1);
        kdata_pa_tmp = epi2D_dis_phs_fw2(x_im,t_es,fm,phs_pa,samp_pa,0);
        
        yy = epi2D_dis_phs_fw2_adj(kdata_ap_tmp,kdata_pa_tmp,t_es,fm,phs_ap, phs_pa, samp_ap,samp_pa);
        
        yy = fft2c(yy);
        
        x_calib = GOP'*(GOP*x_k);
        
        y = [yy(:)+alpha*x_calib(:); lambda*x(1:l)];
    
    end
end