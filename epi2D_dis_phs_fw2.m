function dataout = epi2D_dis_phs_fw2(datain,t_es,fm,phs,sampling_order,ap_tag)
% Forward model for distortion and phase corruption
% assumption : distortion along x is negligible, phase is constant within
% slab
 
[Ny,Nz,Nch] = size(datain);
data_dist = zeros(Ny,Nz,Nch);
dataout = zeros(Ny,Nz,Nch);
nshot = size(sampling_order,3);
% forward - distortion
n_t = Ny;
t_es_e = t_es * size(sampling_order,1) / Ny; % t_es should be the effective echo spacing between adjacent ky lines.
 
for i_t = 1 : n_t
    if ap_tag
        tt = (i_t-0.5) * t_es_e;
    else
        tt = (n_t - i_t + 0.5) * t_es_e;
    end
    ph_accu  = exp(1i * tt * 2*pi * fm);
    ph_accu = repmat(ph_accu,[1,1,Nch]);
    datatmp = fft2c(bsxfun(@times, datain, ph_accu));
    data_dist(i_t,:,:) = datatmp(i_t,:,:);
end
% phase error modulation and sampling
im_data_dist = ifft2c(data_dist); % distorted image
% im_data_dist = datain;

for i_shot = 1 : nshot
    phs_shot = squeeze(phs(:,i_shot));
    
    datatmp = fft2c(bsxfun(@times, im_data_dist, phs_shot));
    ky_curr = sampling_order(:,1,i_shot);
    kz_curr = sampling_order(:,2,i_shot);
    
    indy = ky_curr(find(ky_curr(:)));
    indz = kz_curr(find(kz_curr(:)));
    
    for iii = 1 : length(indy)
        dataout(indy(iii),indz(iii),:)=datatmp(indy(iii),indz(iii),:);
    end
   
end

% adding phase error, and k-space sampling

end
