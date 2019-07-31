
function [DDMobs,L1] = readL1(filename,ddm_index,sample_index)
%read L1 DDM into DDMobs
%ddm_index,sample_index are one based
quality_flags = ncread(filename,'quality_flags');
L1.flags = quality_flags(ddm_index,sample_index);

sp_inc_angle = ncread(filename,'sp_inc_angle');
L1.inc_angle = sp_inc_angle(ddm_index,sample_index);

sp_lat = ncread(filename,'sp_lat');
sp_lon = ncread(filename,'sp_lon');
L1.sp_ll = [sp_lat(ddm_index,sample_index) sp_lon(ddm_index,sample_index)]; % sp lat lon

brcs_ddm_sp_bin_delay_row = ncread(filename,'brcs_ddm_sp_bin_delay_row');
brcs_ddm_sp_bin_dopp_col = ncread(filename,'brcs_ddm_sp_bin_dopp_col');
L1.delay_bin = brcs_ddm_sp_bin_delay_row(ddm_index,sample_index) + 1;  %one based
L1.Doppler_bin = brcs_ddm_sp_bin_dopp_col(ddm_index,sample_index) + 1; %one based

ddm_snr = ncread(filename,'ddm_snr');
L1.SNR = ddm_snr(ddm_index,sample_index);  %-2 to 27

gps_eirp = ncread(filename,'gps_eirp');
L1.eirp = gps_eirp(ddm_index,sample_index);

sp_delay_error = ncread(filename,'sp_delay_error');  %in chips
sp_dopp_error = ncread(filename,'sp_dopp_error');   %in Hz
L1.sp_delay_error = sp_delay_error(ddm_index,sample_index);
L1.sp_dopp_error = sp_dopp_error(ddm_index,sample_index);

Power_analog = ncread(filename,'power_analog');
DDMobs = Power_analog(:,:,ddm_index,sample_index);
DDMobs = DDMobs';%17x11
DDMobs = reshape(DDMobs,[187 1]);
end

