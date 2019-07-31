set(groot,'DefaultFigureColormap',jet);
clear;
addpath('/Users/fax/CYGNSS/VAM/VAM_MATLAB');

%% initialize variables
%choose terms in the cost function
terms.bkg=1;
terms.ddm=1;
terms.L2=0;
terms.div=1;
terms.vor=1;
terms.lap=1;

%scale factor
std_bkg = 10;
scale.bkg=1*1/(std_bkg.^2);
scale.ddm=10;
scale.L2=1/100;
scale.div=1;
scale.vor=1;
scale.lap=2;

%choose DDM bins to use: 'specular','full','rectangle','threshold','amb_free'
DDMtype = 'threshold';

%DDM covariance matrix: 'constant', 'scale', 'numerical'
load iR;
DDMcov.type = 'scale';
DDMcov.constant = 1e35; %for 'constant'
DDMcov.iR = iR;
DDMcov.scale = 0.10; %for 'scale' type

%parameters to stop iteration
Iter.Tol=1e-3;
Iter.MaxIter=10;
Iter.MaxFun=15;

%Quality Control
QC.inc_angle_max = 60;
QC.SNR_min = 3;
QC.ws_min = 3;
QC.ws_max = 35;
QC.diff_max = 1;  %0.5
QC.correlation_min = 0.9; %0.95 

%% prepare the initial background
% BKG
% First choose a fixed background area (can be global), then choose the integration area
% The integration area is around the specular point
% For now, all have resolution of 0.125 degree

%Big map: Nlat, Nlon, Npts, LON_vec, LAT_vec, U, V, UV_bkg, UV_ana
resolution = 0.125;
LON_min = 301.1;
LAT_min = 13;
Nlon = 49; %grid points of the large map
Nlat = 49;
Npts=Nlat*Nlat; %number of points for the large map; Nstate = Npts*2
wind = '../../Data/Irma2017/irma11l.2017090418.hwrfprs.synoptic.0p125.f005.uv.nc';
LON_vec=ncread(wind,'longitude'); 
LAT_vec=ncread(wind,'latitude'); 
U=ncread(wind,'UGRD_10maboveground');
V=ncread(wind,'VGRD_10maboveground');
iLON = find(abs(LON_vec-LON_min)<0.01):find(abs(LON_vec-LON_min)<0.01)+Nlon-1;
iLAT = find(abs(LAT_vec-LAT_min)<0.01):find(abs(LAT_vec-LAT_min)<0.01)+Nlat-1;
BKG.LON_vec=LON_vec(iLON);
BKG.LAT_vec=LAT_vec(iLAT);
BKG.U=U(iLON,iLAT); %49x49
BKG.V=V(iLON,iLAT); %49x49

% WS_BKG = sqrt(BKG.U.^2+BKG.V.^2);
% figure;imagesc(BKG.LON_vec,BKG.LAT_vec,WS_BKG');colorbar;set(gca,'YDir','normal');

%% CYGNSS
CYGNSS.filename = '/users/fax/CYGNSS/Data/Irma2017/cyg04.ddmi.s20170904-000000-e20170904-235959.l1.power-brcs.a21.d21.nc';
CYGNSS.ddm_index = 1;%one based 1-4

%% start of assimilation
int_size = 17; %intergration size; 17x0.125 = 2.125deg; must be larger than 120km x 120km; must be odd number
temp_path = './temp/'; %directory for config file and other temp files; may need to mkdir for it
ANA = BKG;
for index = [81045 81096] %81040-81100
    CYGNSS.index = index;
    ANA = vam_main(terms,scale,DDMtype,DDMcov,Iter,QC,ANA,CYGNSS,temp_path,int_size);
end

%% plot map
ANA.WS = sqrt(ANA.U.^2+ANA.V.^2);
BKG.WS = sqrt(BKG.U.^2+BKG.V.^2);

figure('units','normalized','outerposition',[0.02 0.6 0.95 0.4]); %left bottom width height
subplot(1,3,1);
imagesc(BKG.LON_vec,BKG.LAT_vec,BKG.WS');colorbar;
xlabel('Longitude','FontSize',14);ylabel('Latitude','FontSize',14);set(gca,'YDir','normal');set(gca,'Fontsize',16)
subplot(1,3,2);
imagesc(ANA.LON_vec,ANA.LAT_vec,ANA.WS');colorbar;
xlabel('Longitude','FontSize',14);ylabel('Latitude','FontSize',14);set(gca,'YDir','normal');set(gca,'Fontsize',16)
subplot(1,3,3);
imagesc(ANA.LON_vec,ANA.LAT_vec,ANA.WS'-BKG.WS');colorbar;
xlabel('Longitude','FontSize',14);ylabel('Latitude','FontSize',14);set(gca,'YDir','normal');set(gca,'Fontsize',16)


