function [J,g] = costFun(size,uv_ana,uv_bkg,terms,scale,DDM,L2,temp_path,fm_path)
%VAM cost function including background terms and observation terms

%size is the number of points in one side of the region (might use global variable)
%uv_ana & uv_bkg are vector 1xn
%DDMobs is a 187x1 vector (can be reshape into 17x11)
%bin_index is the index vector of effective bins
%iR is 187x187 inverse covariance matrix of DDM

%J = scale_bkg*sum(uv_ana-uvbkg)^2 + scale_ddm*(DDMobs-DDMfm)^2

resolution=0.125;

n=size*size*2; 
Jb=0;Jddm=0;Jdiv=0;Jvor=0;Jl2=0;Jlap=0;
gb=zeros(n,1);
gddm=zeros(n,1);
gl2=zeros(n,1);
gdiv=zeros(n,1);
gvor=zeros(n,1);
glap=zeros(n,1);

%% Background terms
if (terms.bkg==1)
    Jb=scale.bkg*sum((uv_ana-uv_bkg).^2);
    gb=2*scale.bkg*(uv_ana-uv_bkg)';  %nx1 gradient for background terms
end

%% DDM terms
if (terms.ddm==1)
    %preapare input wind for forward model
    fid1=fopen([temp_path,'uv_input.dat'],'w');
    fwrite(fid1,uv_ana,'double');
    fclose(fid1);

    %run forward model
    [status,cmdout]=system([fm_path,' ',temp_path,'config.txt']); %status=0 means the command completed successfully
    if(status~=0)
        disp('The forward model does not run successfully with cmdout:');
        disp(cmdout);
    end
    %disp(cmdout);
    
    %read DDMfm
    fid2 = fopen([temp_path,'DDMfm.dat']);
    DDMfm= fread(fid2,'double'); %read DDMfm 187x1
    fclose(fid2);
    
    bin_index=DDM.bin_index;
    Nbin=length(bin_index);
    DDMfm=DDMfm(bin_index);
    DDMobs=DDM.power(bin_index);
    iR=DDM.iR(bin_index,bin_index);
    Jddm=scale.ddm*(DDMfm-DDMobs)'*iR*(DDMfm-DDMobs);

    %gradient
    fid1 = fopen([temp_path,'Jacobian.dat']);
    data = fread(fid1,'double')';
    fclose(fid1);
    numDDMbins=data(1);
    numSurfacePts=data(2);
    H=reshape(data(3:end),[numDDMbins,numSurfacePts]); %read H: 187x128
    fid1 = fopen([temp_path,'indexLL.dat']);
    indexLL = fread(fid1,'int');  %index of the 128 points in the grid n/2; zeros based
    indexLL=indexLL+1; %change to 1 based
    fclose(fid1);
    
    H=H(bin_index,:); %Nbinx128
    u=uv_ana(1:n/2)'; %n/2x1
    v=uv_ana(n/2+1:n)'; %n/2x1
    ws=sqrt(u.^2+v.^2); %n/2x1
    gu=zeros(n/2,Nbin);
    gv=zeros(n/2,Nbin);
    gu(indexLL,:)=H';
    gv(indexLL,:)=H';
    for i=1:Nbin
        gu(:,i)=gu(:,i).*(u./ws); %chain rule: from wind speed to u,v
        gv(:,i)=gv(:,i).*(v./ws);
    end
    H_uv=[gu;gv];%nx187
    gddm=scale.ddm*2*H_uv*iR*(DDMfm-DDMobs);
end

%% L2 wind speed term
if (terms.L2==1)
    % L2.ws L2.lat L2.lon

    %interpolate uv_ana to the wind speed observation location
    u=uv_ana(1:n/2)';
    v=uv_ana(n/2+1:n)'; %n/2x1
    ws=sqrt(u.^2+v.^2); %n/2x1
    ws_2d=reshape(ws,[size size]);
    ws_ana=interp2(L2.lat_vec,L2.lon_vec,ws_2d,L2.lat,L2.lon,'linear');
    Jl2=scale.L2*(ws_ana-L2.ws)^2;
    
    %gradient
    %ws_ana=ws(bi_index(1))*bi_weight(1)+...+ws(bi_index(4))*bi_weight(4)
    %ws=sqrt(u^2+v^2) .  dws=u/w*du+v/w*dv
    [bi_index, bi_weight] = bilinear_interp(lon_vec,lat_vec,L2.lon,L2.lat,resolution);

    gu=zeros(n/2,1);
    gv=zeros(n/2,1);
    gu(bi_index)=bi_weight;
    gv(bi_index)=bi_weight;
    gu=gu.*(u./ws);
    gv=gv.*(v./ws);
    gl2=2*scale.L2*(ws_ana-L2.ws)*[gu;gv];
end

%% Divergence term
if (terms.div==1)
    [Jdiv,gdiv]=Jdivergence(size,uv_ana,uv_bkg,scale.div);
%     %finite difference to compute gradient
%     for i=1:n
%         uv_temp=uv_ana;
%         uv_temp(i)=uv_ana(i)+1e-10;
%         Jdiv2=Jdivergence(size,uv_temp,uv_bkg,scale.div);
%         gdiv(i)=(Jdiv2-Jdiv)/1e-10;
%     end
end

%% Vorticity term
if (terms.div==1)
    [Jvor,gvor]=Jvorticity(size,uv_ana,uv_bkg,scale.vor);
end

%% Laplacian term
if (terms.lap==1)
    [Jlap,glap]=Jlaplacian(size,uv_ana,uv_bkg,scale.lap);
end

%%
J=Jb+Jddm+Jl2+Jdiv+Jvor+Jlap;
g=gb+gddm+gl2+gdiv+gvor+glap;
disp(['J = ',num2str(J),', Jb = ',num2str(Jb),', Jddm = ',num2str(Jddm), ', Jl2 = ',num2str(Jl2),...
    ', Jdiv = ', num2str(Jdiv), ', Jvor = ', num2str(Jvor), ', Jlap = ', num2str(Jlap)]);
end

