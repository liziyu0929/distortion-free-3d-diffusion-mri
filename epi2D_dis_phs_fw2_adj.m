function dataout = epi2D_dis_phs_fw2_adj(datain_ap,datain_pa,t_es,fm,phs_ap, phs_pa, samp_ap,samp_pa)
% adjoint of forward model for distortion and phase corruption
% assumption : distortion along x is negligible, phase is constant within
% slab

[Ny,Nz,Nch] = size(datain_ap);
im_dist_ap = zeros(Ny,Nz,Nch);
im_dist_pa = zeros(Ny,Nz,Nch);
dataout = zeros(Ny,Nz,Nch);

n_shot = size(samp_ap,3);

for i_shot = 1 : n_shot
    
    phs_ap_shot = squeeze(phs_ap(:,i_shot));
    phs_pa_shot = squeeze(phs_pa(:,i_shot));
    
    data_ap_tmp = zeros(Ny,Nz,Nch);
    data_pa_tmp = zeros(Ny,Nz,Nch);
    
    ky_curr = samp_ap(:,1,i_shot);
    kz_curr = samp_ap(:,2,i_shot);
    indy = ky_curr(find(ky_curr(:)));
    indz = kz_curr(find(kz_curr(:)));
    
    for iii = 1 : length(indy)
        data_ap_tmp(indy(iii),indz(iii),:)=datain_ap(indy(iii),indz(iii),:);
    end
    
    im_ap_tmp = bsxfun(@times, ifft2c(data_ap_tmp), conj(phs_ap_shot));
    im_dist_ap = im_dist_ap + im_ap_tmp;
    
    ky_curr = samp_pa(:,1,i_shot);
    kz_curr = samp_pa(:,2,i_shot);
    indy = ky_curr(find(ky_curr(:)));
    indz = kz_curr(find(kz_curr(:)));
    
    for iii = 1 : length(indy)
        data_pa_tmp(indy(iii),indz(iii),:)=datain_pa(indy(iii),indz(iii),:);
    end
    
    im_pa_tmp = bsxfun(@times, ifft2c(data_pa_tmp), conj(phs_pa_shot));
    im_dist_pa = im_dist_pa + im_pa_tmp;

end

% dataout = im_dist;

data_dist_ap = fft2c(im_dist_ap);
data_dist_pa = fft2c(im_dist_pa);

n_t = Ny;

t_es_e = t_es * size(samp_ap,1) / Ny; % t_es should be the effective echo spacing between adjacent ky lines.


for i_t = 1 : n_t
    
    data_ap_tmp = zeros(Ny,Nz,Nch);
    data_pa_tmp = zeros(Ny,Nz,Nch);
    
    tt_ap = (i_t-0.5) * t_es_e;
    tt_pa = (n_t - i_t + 0.5) * t_es_e;
    
    ph_accu_ap  = exp(1i * tt_ap * 2*pi * fm);
    ph_accu_pa  = exp(1i * tt_pa * 2*pi * fm);
    
    ph_accu_ap = repmat(ph_accu_ap,[1,1,Nch]);
    ph_accu_pa = repmat(ph_accu_pa,[1,1,Nch]);
    
    data_ap_tmp(i_t,:,:) = data_dist_ap(i_t,:,:);
    data_pa_tmp(i_t,:,:) = data_dist_pa(i_t,:,:);
    
    im_ap_tmp = bsxfun(@times, ifft2c(data_ap_tmp), conj(ph_accu_ap));
    im_pa_tmp = bsxfun(@times, ifft2c(data_pa_tmp), conj(ph_accu_pa));
        
    dataout = dataout + im_ap_tmp + im_pa_tmp;
        
end


end




