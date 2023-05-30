
clc, clear,close all

Rz = 2; nIterCG = 30; lambda = 2e-4; alpha = 20;
L1 = 0; % 0 for recon w/o L1 constraint, 1 for recon w/ L1
beta = 0.4; nIterL1 = 12; nIterSplit = 4; lambda2 = 5e-7;
pf = 1;

% lambda: Tiknov regularization weights; alpha: SPIRiT constraint weights
% nIterCG: iter number for joint recon; nIterL1: iter number for L1
% refinement; nIterSplit: number of split iter for L1; beta: split weight
% for L1; lambda2: L1 weight


addpath(genpath('utils'));

nx = 180; ny = 180; nz = 25; nch = 8;
t_es = 0.76e-3;

% calibration data
load data/gre_calib.mat
data_calib = fft3c(im_calib);

if Rz == 2
    if pf == 1
        ap_index = [1:2:9, 11:15];
        pa_index = [11:15, 16:2:24];
    else
        ap_index = 1:2:24;
        pa_index = 2:2:24;
    end
else
    ap_index = 1:24;
    pa_index = 1:24;
end

% AP data
load('data/raw_dwi_ap_1.22_25sli.mat');
load('data/nav_phs_ap.mat');
kAP = permute(squeeze(k_data_3DMS_AP),[1,2,4,3]);
phs_ap = pherrall_AP(:,:,ap_index);
samp_ap = samp_all_AP(:,:,ap_index);

% PA data
load('data/raw_dwi_pa_1.22_25sli.mat');
load('data/nav_phs_pa.mat');
kPA = permute(squeeze(k_data_3DMS_PA),[1,2,4,3]);
phs_pa = pherrall_PA(:,:,pa_index);
samp_pa = samp_all_PA(:,:,pa_index);

% field map
load('data/field_map.mat')
fm = fm_caipi;

% reference data
load('data/dwi_ref_1.22_25sli.mat');
imref = permute(abs(imref_joint),[2,1,3]);
imref = imrecon_out(end:-1:1,:,:);
figure,montage(imrecon_out(:,:,4:end-2),'Size',[4,5],'DisplayRange',[0 max(imrecon_out(:))*0.3]);

% coil compression
kdata_ap = coil_compression(kAP,data_calib,nch);
kdata_pa = coil_compression(kPA,data_calib,nch);
data_calib = coil_compression(data_calib,data_calib,nch);

mask_yz_ap = sum(samp_ap, 3);
mask_ap = repmat(reshape(mask_yz_ap, [1 size(mask_yz_ap)]), [nx, 1, 1]);
kdata_ap_ds = bsxfun(@times, kdata_ap, mask_ap);

mask_yz_pa = sum(samp_pa, 3);
mask_pa = repmat(reshape(mask_yz_pa, [1 size(mask_yz_pa)]), [nx, 1, 1]);
kdata_pa_ds = bsxfun(@times, kdata_pa, mask_pa);

samp_ap = shot2samp(samp_ap, 2.5, 3);
samp_pa = shot2samp(samp_pa, 2.5, 3);

kdata_ap_ds_1 = ifft1c(kdata_ap_ds,1);
kdata_pa_ds_1 = ifft1c(kdata_pa_ds,1);
data_calib_fro = ifft1c(data_calib,1);


% SPIRiT
kSize = [5,5];
CalibTyk = 0.1;
sz=[ny,nz];
imrecon = zeros(nx,ny,nz,nch);

parfor ix = 1:nx
    disp(['ix:' num2str(ix)])
    
    kdata_ap_x =  squeeze(kdata_ap_ds_1(ix,:,:,:));
    kdata_pa_x =  squeeze(kdata_pa_ds_1(ix,:,:,:));
    phs_ap_x = conj(squeeze(phs_ap(ix,:,:)));
    phs_pa_x = conj(squeeze(phs_pa(ix,:,:)));
    fm_x = squeeze(fm(ix,:,:));

    DATA_x_calib = squeeze(data_calib_fro(ix,:,:,:));
    kernel = calibSPIRiT(DATA_x_calib, kSize, nch, CalibTyk);
    GOP = SPIRiT(kernel, 'fft', [ny, nz]);

    [resSPIRiT] = pcgSPIRiT_distortionPhaseCorrection2(...
        kdata_ap_x, kdata_pa_x,fm_x,t_es, phs_ap_x, phs_pa_x,...
        samp_ap, samp_pa, GOP, nIterCG, sz, nch, lambda, alpha);

    if L1 == 1
        XOP = Wavelet('Daubechies_TI',4,6);
        [resSPIRiT] = pcgSPIRiTL1_distortionPhaseCorrection(...
            kdata_ap_x, kdata_pa_x,fm_x,t_es, phs_ap_x, phs_pa_x,...
            samp_ap, samp_pa, GOP, XOP, nIterL1, nIterSplit, sz, nch, lambda2,...
            alpha, beta, resSPIRiT);
    end

    imrecon(ix,:,:,:) = squeeze(ifft2c(resSPIRiT));
    
end

imrecon_out = permute(sos(imrecon),[2,1,3]);
imrecon_out = imrecon_out(end:-1:1,:,:);
figure,montage(imrecon_out(:,:,4:end-2),'Size',[4,5],'DisplayRange',[0 max(imrecon_out(:))*0.3]);