function [res,RESVEC] = pcgSPIRiT_phaseCorrection(kdata, phs, sampling_order, GOP, nIter, sz, nch,z_start,lambda,alpha)

y = epi2D_phs_fw_adj(kdata,phs,sampling_order,z_start);
y = fft2c(y);

N = prod(sz)*nch;
z = zeros([sz nch]); % Tikhonov regularization term
b = [y(:); z(:)];

[tmpres,FLAG,RELRES,ITER,RESVEC] = pcg(@(x)afun(x,phs,sampling_order,GOP,sz,nch,z_start,lambda,alpha),b,1e-6,nIter);


res = reshape(tmpres(1:N),[sz nch]);

    function y = afun(x,phs,sampling_order,GOP,sz,nch,z_start,lambda,alpha)
        
        l = prod(sz)*nch;
        x_k_in = reshape(x(1:l), [sz nch]);
        x_im = ifft2c(x_k_in);
        x_k = epi2D_phs_fw(x_im, phs, sampling_order,z_start); % DP
        x_im = epi2D_phs_fw_adj(x_k,phs,sampling_order,z_start); % (DP)^H
        x_k = fft2c(x_im);
        x_calib = GOP'*(GOP*x_k_in);
        y = [x_k(:)+alpha*x_calib(:); lambda*x(1:l)];
    end
end