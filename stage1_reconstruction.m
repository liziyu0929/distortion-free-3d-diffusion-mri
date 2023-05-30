clc, clear, close all

Rz = 2; pf = 1; nIterCG = 30; mu = 2e-4; lambda1 = 1;
L1 = 0; % 0 for recon w/o L1 constraint, 1 for recon w/ L1
beta = 0.4; nIterL1 = 12; nIterSplit = 4; lambda2 = 5e-7;
bu = 1;

% mu: Tiknov regularization weights; alpha: SPIRiT constraint weights
% nIterCG: iter number for joint recon; nIterL1: iter number for L1
% refinement; nIterSplit: number of split iter for L1; beta: split weight
% for L1; lambda2: L1 weight; bu: 1-blip-up; 0-blip-down

%% data and parameters

addpath(genpath('utils'));

% calibration data
load data/gre_calib.mat
data_calib = fft3c(im_calib);

% AP data
if Rz == 2
    if pf == 1
        if bu == 1
            index = [1:2:9, 11:15];
        else
            index = [11:15, 16:2:24];
        end
    else
        index = 1:2:24;
    end
else
    index = 1:24;
end

if bu == 1
    load('data/raw_dwi_ap_1.22_25sli.mat');
    load('data/nav_phs_ap.mat');
    kdata = permute(squeeze(k_data_3DMS_AP),[1,2,4,3]);
    phs = pherrall_AP(:,:,index);
    samp = samp_all_AP(:,:,index);
else
    load('data/raw_dwi_pa_1.22_25sli.mat');
    load('data/nav_phs_pa.mat');
    kdata = permute(squeeze(k_data_3DMS_PA),[1,2,4,3]);
    phs = pherrall_PA(:,:,index);
    samp = samp_all_PA(:,:,index);
end

[nx, ny, nz, nch] = size(kdata);

% coil compression
nch = 8;
data3d = coil_compression(kdata, data_calib, nch);
data_calib = coil_compression(data_calib, data_calib, nch);

mask_yz = sum(samp, 3);
mask = repmat(reshape(mask_yz, [1 size(mask_yz)]), [nx, 1, 1]);
data_ds = bsxfun(@times, data3d, mask);

if Rz == 2 && pf == 1
    if bu == 1
        data_ds1 = data_ds(:,:,1:17,:);
        z_start = 1;
    else
        data_ds1 = data_ds(:,:,9:end,:);
        z_start = 9;
    end
else
    data_ds1 = data_ds;
end

sampling_order = shot2samp(samp, 2, 3);

% SPIRiT
kSize = [5,5];
CalibTyk = 0.1;

kdata_ds_1 = ifft1c(data_ds1,1);
data_calib_fro = ifft1c(data_calib,1);
nz_eff = size(kdata_ds_1, 3);
sz = [ny,nz_eff];

imrecon = zeros(nx,ny,nz,nch);

parfor ix = 1:nx
    disp(['ix:' num2str(ix)])
    kdata_x = squeeze(kdata_ds_1(ix,:,:,:));
    data_phs_x = conj(squeeze(phs(ix,:,:)));
    DATA_x_calib = squeeze(data_calib_fro(ix,:,:,:));
    kernel = calibSPIRiT(DATA_x_calib, kSize, nch, CalibTyk);
    GOP = SPIRiT(kernel, 'fft', [ny, nz_eff]);
    [resSPIRiT] = pcgSPIRiT_phaseCorrection(...
        kdata_x, data_phs_x, sampling_order, GOP, nIterCG, sz, nch, z_start, mu, lambda1);
    
    if L1 == 1
        XOP = Wavelet('Daubechies_TI',4,6);
        [resSPIRiT] = pcgSPIRiTL1_phaseCorrection(...
            kdata_x, data_phs_x, sampling_order, GOP, XOP,...
            nIterL1, nIterSplit, sz, nch, z_start, lambda2,...
            lambda1, beta, resSPIRiT);
    end
    
    if Rz == 2 && pf == 1
        if bu == 1
            k_out = cat(2,squeeze(resSPIRiT),squeeze(data_ds(ix,:,18:end,:)));
        else
            k_out = cat(2,squeeze(data_ds(ix,:,1:8,:)),squeeze(resSPIRiT));
        end
    else
        k_out = squeeze(resSPIRiT);
    end
    
    imrecon(ix,:,:,:) = squeeze(ifft2c(k_out));
end

imrecon_out = permute(sos(imrecon),[2,1,3]);
imrecon_out = imrecon_out(end:-1:1,:,:);
figure,montage(imrecon_out(:,:,4:end-2),'Size',[4,5],'DisplayRange',[0 max(imrecon_out(:))*0.3]);