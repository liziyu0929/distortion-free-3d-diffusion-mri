function samp_order = shot2samp(samp_shot, Rz, Ry)
[ny, nz, nshot] = size(samp_shot);
% samp_order = zeros(ceil(ny/Ry),2,1);

for ii = 1 : nshot
    idx_samp = find(samp_shot(:,:,ii)==1);
    [ky, kz] = ind2sub(size(samp_shot(:,:,ii)), idx_samp);
    [ky_s, idx] = sort(ky);
    kz_s = kz(idx);
    samp_order(:, 1, ii) = ky_s; % ky sampling
    samp_order(:, 2, ii) = kz_s; % kz sampling
end
end
    