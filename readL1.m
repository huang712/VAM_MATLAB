function L1 = readL1(filename,ddm_index,sample_index)
% Read CYGNSS L1 data
% ddm_index,sample_index are one based

L1.filename = filename;
L1.ddm_index = ddm_index;
L1.index = sample_index;

L1.sample_time=ncread(filename,'ddm_timestamp_utc',sample_index,1);
L1.flags=ncread(filename,'quality_flags',[ddm_index sample_index],[1 1]);
L1.inc_angle=ncread(filename,'sp_inc_angle',[ddm_index sample_index],[1 1]);
L1.sp_ll(1)=ncread(filename,'sp_lat',[ddm_index sample_index],[1 1]); % lat
L1.sp_ll(2)=ncread(filename,'sp_lon',[ddm_index sample_index],[1 1]); % lon
L1.delay_bin=ncread(filename,'brcs_ddm_sp_bin_delay_row',[ddm_index sample_index],[1 1])+1; % one based
L1.Doppler_bin=ncread(filename,'brcs_ddm_sp_bin_dopp_col',[ddm_index sample_index],[1 1])+1; % one based
L1.SNR=ncread(filename,'ddm_snr',[ddm_index sample_index],[1 1]); % -2 to 27
L1.eirp=ncread(filename,'gps_eirp',[ddm_index sample_index],[1 1]); % Watt
L1.sp_delay_error=ncread(filename,'sp_delay_error',[ddm_index sample_index],[1 1]); % in chips
L1.sp_dopp_error=ncread(filename,'sp_dopp_error',[ddm_index sample_index],[1 1]); % in Hz
L1.scNum=ncread(filename,'spacecraft_num');
L1.prn_code=ncread(filename,'prn_code',[ddm_index sample_index],[1 1]);
L1.ddm_ant=ncread(filename,'ddm_ant',[ddm_index sample_index],[1 1]);
L1.sp_rx_gain=ncread(filename,'sp_rx_gain',[ddm_index sample_index],[1 1]);

DDMobs = ncread(filename,'power_analog',[1 1 ddm_index sample_index],[11 17 1 1]);
DDMobs=DDMobs'; % 17x11
L1.DDMobs = reshape(DDMobs,[187 1]);

end

