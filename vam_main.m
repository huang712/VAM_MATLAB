 %Main function of VAM

function [ANA,flag] = vam_main(terms,scale,DDMtype,DDMcov,Iter,QC,BKG,L1,temp_path,fm_path,int_size)
ANA = BKG;
flag = 1; %when pass QC, flag = 0

%First Quality control 
%----------------------------------------------------------------------
if (qfCYG(L1.flags)~=0) 
    disp('Poor overall quality');
    return;
end

if (L1.inc_angle > QC.inc_angle_max) 
    disp(['Large incidence angle: ',num2str(L1.inc_angle)]);
    return;
end

if (L1.SNR < QC.SNR_min) 
    disp(['Low SNR: SNR = ',num2str(L1.SNR)]);
    return;
end

%check if specular point(integration area) is inside the map
splat_index=find(abs(BKG.LAT_vec-L1.sp_ll(1))==min(abs(BKG.LAT_vec-L1.sp_ll(1))));
splon_index=find(abs(BKG.LON_vec-L1.sp_ll(2))==min(abs(BKG.LON_vec-L1.sp_ll(2))));
splat_index=splat_index(1);splon_index=splon_index(1);
k=(int_size-1)/2;
if(splat_index-k<1 || splat_index+k>length(BKG.LAT_vec) || splon_index-k<1 || splon_index+k>length(BKG.LON_vec))
    disp('Specular point is outside the map');return;
end

%QC on wind speed
%find wind speed of background at specular point
WS = sqrt(BKG.U.^2+BKG.V.^2);
sp_ws=WS(splon_index,splat_index);
if(sp_ws < QC.ws_min || sp_ws > QC.ws_max)
    disp(['wind speed = ',num2str(sp_ws),' does not pass wind speed quality control']);
    return;
end
%----------------------------------------------------------------------

%Use the specular point to compute the integration area 
%Intergration area: nlat, nlon, npts, lon_vec, lat_vec, u, v, uv_bkg, uv_ana
%----------------------------------------------------------------------
resolution = 0.125;
nlat = int_size; nlon = int_size;
npts = nlat*nlon; %number of state in the integration area
lat_vec = BKG.LAT_vec(splat_index-k:splat_index+k);
lon_vec = BKG.LON_vec(splon_index-k:splon_index+k);
u_bkg = BKG.U(splon_index-k:splon_index+k,splat_index-k:splat_index+k);
v_bkg = BKG.V(splon_index-k:splon_index+k,splat_index-k:splat_index+k);
uv_bkg = zeros(1,npts*2);
uv_bkg(1:npts)=u_bkg; uv_bkg(npts+1:npts*2)=v_bkg;
uv_ana = uv_bkg; 

% figure;imagesc(BKG.LON_vec(splon_index-k:splon_index+k),BKG.LAT_vec(splat_index-k:splat_index+k),WS(splon_index-k:splon_index+k,splat_index-k:splat_index+k)');set(gca,'YDir','normal');
% ws_bkg = sqrt(u_bkg.^2+v_bkg.^2);
% figure;imagesc(lon_vec,lat_vec,ws_bkg');set(gca,'YDir','normal');hold on
% plot(L1.sp_ll(2),L1.sp_ll(1),'r*');
%----------------------------------------------------------------------

%write index to config file for Forward model to read
fid=fopen([temp_path,'config.txt'],'w');
fprintf(fid,'%s\n',[temp_path,'uv_input.dat']); %wind field data
fprintf(fid,'%s\n',L1.filename);
fprintf(fid,'%s\n',temp_path);  %path to save DDMfm, Jacobian and indexLL
fprintf(fid,'ddmIndex    = %d\n',L1.ddm_index-1);  %zero based
fprintf(fid,'sampleIndex = %d\n',L1.index-1);  %zero based
fprintf(fid,'numPtsLon   = %d\n',nlon);  
fprintf(fid,'numPtsLat   = %d\n',nlat);
fprintf(fid,'lon_min_deg = %.3f\n',BKG.LON_vec(splon_index-k));
fprintf(fid,'lat_min_deg = %.3f\n',BKG.LAT_vec(splat_index-k));
fprintf(fid,'resolution  = %.3f\n',resolution);
fprintf(fid,'JacobOnoff  = %d\n',1);
fprintf(fid,'thermalNoiseOnOff = %d\n',0);
fclose(fid);

%write DDMobs to file
DDMobs = L1.DDMobs; %187x1
fid=fopen([temp_path,'DDMobs.dat'],'w');
fwrite(fid,L1.DDMobs,'double');
fclose(fid);

%choose DDM bins
if (strcmp(DDMtype,'specular'))
    bin_index = (round(L1.Doppler_bin)-1)*17 + round(L1.delay_bin);
elseif (strcmp(DDMtype,'full'))
    bin_index = 1:187;
elseif (strcmp(DDMtype,'rectangle')) %11x5 bins
    a = reshape(1:187,[17 11]);
    bin_index = a(round(L1.delay_bin):17,4:8);   
    bin_index = reshape(bin_index,[1 length(bin_index(:))]);
elseif (strcmp(DDMtype,'threshold')) %use 1/10 peak power as threshold
    a = reshape(1:187,[17 11]);
    info_bins=reshape(a(4:17,3:9),[1 14*7]); %discard first three rows
    bin_index = find(DDMobs>max(DDMobs)/10);
    bin_index = intersect(bin_index,info_bins);
elseif (strcmp(DDMtype,'amb_free'))
    bin_index = [48 62 76 92 110 130 150];
    disp('Ambiguity free line mode is not working now');
    %return;
end

%DDM covariance
if (strcmp(DDMcov.type,'constant'))
    DDM.iR = DDMcov.constant*eye(187);
elseif (strcmp(DDMcov.type,'scale'))
    DDM.iR = diag(1./((DDMobs*DDMcov.scale).^2));
    DDM.iR = DDM.iR(bin_index,bin_index);
elseif (strcmp(DDMcov.type,'numerical'))
    DDM.iR = DDMcov.iR;
elseif (strcmp(DDMcov.type,'model'))
    DDM.iR = computeCov(DDMobs,bin_index,L1.delay_bin,L1.Doppler_bin,sp_ws,DDMcov.method,DDMcov.k_max);
end

%load DDM observation
DDM.power=DDMobs;
DDM.bin_index=bin_index;

%CYGNSS L2 wind speed observation
L2.ws=17.9;
L2.lat=17.24;
L2.lon=304.87;
L2.lon_vec = lon_vec;
L2.lat_vec = lat_vec;

%compute initial cost function
%disp('Compute initial cost function')

[J0,g0]=costFun(int_size,uv_ana,uv_bkg,terms,scale,DDM,L2,temp_path,fm_path);
movefile([temp_path,'DDMfm.dat'],[temp_path,'DDMfm0.dat'])  %simulated DDM from background

%Second quality control: compare DDMobs and DDMfm0 187*1
%----------------------------------------------------------------------
%1. absolute average relative power difference  2. correlation coefficient
fid=fopen([temp_path,'DDMfm0.dat']);
DDMfm0 = fread(fid,'double'); %read DDMfm 187x1
fclose(fid);

%relative power difference
diff_DDM1 = (DDMobs(bin_index)-DDMfm0(bin_index))./DDMobs(bin_index);
diff_DDM2 = (DDMfm0(bin_index)-DDMobs(bin_index))./DDMfm0(bin_index);
diff1 = mean(abs(diff_DDM1)); 
diff2 = mean(abs(diff_DDM2)); 
diff = max(diff1,diff2);

%correlation coefficient
correlation = corr2(DDMobs(bin_index),DDMfm0(bin_index)); 
disp(['Relative power difference = ',num2str(diff),'; correlation = ',num2str(correlation)]);
if (diff > QC.diff_max || correlation < QC.correlation_min)
    disp('Does not pass background check: ')
    return;
end
%figure;subplot(1,2,1);imagesc(reshape(DDMobs,[17 11]));colorbar;
%subplot(1,2,2);imagesc(reshape(DDMfm0,[17 11]));colorbar;
%----------------------------------------------------------------------

options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton',...
    'GradObj','on','MaxIterations',Iter.MaxIter,'MaxFunctionEvaluations',Iter.MaxFun,'TolX',Iter.Tol,'TolFun',Iter.Tol); %BFGS with gradient

disp('start minimization:')
tic
[uv_ana,Costval,exitflag,output]=fminunc(@(uv_ana) costFun(int_size,uv_ana,uv_bkg,terms,scale,DDM,L2,temp_path,fm_path),...
    uv_ana,options); %minimization; gives the final uv_ana
toc
disp(output);
disp(['Costval = ',num2str(Costval)]);

%put uv_ana back into the big map UV_ana
u_ana = reshape(uv_ana(1:npts),[nlon,nlat]);
v_ana = reshape(uv_ana(npts+1:npts*2),[nlon,nlat]);
ANA.U(splon_index-k:splon_index+k,splat_index-k:splat_index+k) = u_ana;
ANA.V(splon_index-k:splon_index+k,splat_index-k:splat_index+k) = v_ana;

%gross specular wind speed on Background and Analysis -- this is not the final one!
WS_ana = sqrt(ANA.U.^2+ANA.V.^2);
ws_sp_ana = WS_ana(splon_index,splat_index); %wind speed of analysis at specular point          
disp(['Background wind speed = ',num2str(sp_ws),'; Analysis wind speed = ',num2str(ws_sp_ana)]);

flag = 0;

end

