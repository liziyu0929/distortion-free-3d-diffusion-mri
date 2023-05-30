function data_cc = coil_compression(datain, data_calib, nch_out)

% coil compression using ESPIRiT package
disp('coil compression')
tic
gccmtx = calcGCCMtx(data_calib,1,5);
gccmtx_aligned = alignCCMtx(gccmtx(:,1:nch_out,:));
data_cc = CC(datain,gccmtx_aligned,1);
toc

end
