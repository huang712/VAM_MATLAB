function iR = computeCov_original(DDMobs,bin_index,sp_delay_bin,sp_doppler_bin,sp_ws,method,k_max)
% Old version: using Covariance_parameter.mat (from ECMWF operational forecast)
% compute DDMcovariance matrix using the empirical model
% iR: inverse of the DDM covariance matrix
% DDMobs in vector
% bin_index: index of 187[17x11] of effective bins
% sp_ws: wind speed

load('Covariance_parameter.mat','sigma2_n','A','correlationA','A1','B1','C1');
%delay=-1:10 doppler=-3:3 in parameters
%A 12x7 
%A1,B1,C1: 84x84

N=length(bin_index); %number of effective bins
R=zeros(N,N);

%% find delay,doppler of the effective bins and index (X,Y) in the matrix 12x7
delay=zeros(1,N);
doppler=zeros(1,N);
X=zeros(1,N);
Y=zeros(1,N);
for i=1:N
    [delay(i),doppler(i)]=index2delayDoppler(bin_index(i));
    x=delay(i)-sp_delay_bin+2; %1:12
    y=doppler(i)-sp_doppler_bin+4; %1:7
    if (y<1) 
        y=1;
    elseif (y>7) 
        y=7;
    end
    if (x<1) 
        x=1;
    elseif (x>12) 
        x=12;
    end   
    X(i)=x;
    Y(i)=y;
end

% figure;imagesc(reshape(DDMobs,[17,11]));hold on;
% plot(doppler,delay,'ro')
% plot(sp_doppler_bin,sp_delay_bin,'y*')
% rectangle('Position',[sp_doppler_bin-3.5 sp_delay_bin-1.5 7 12],'EdgeColor','r','LineWidth',2)
% xlabel('delay bin');ylabel('doppler bin');
% set(gca,'fontsize',16);

%% Compute diagonal values
% sigma2=(aY)^2+sigma2_n
a=zeros(1,N);

%find a in A(12,7) using delay,doppler and sp index
for i=1:N
    %find index of [12,8] with respect to sp_delay_bin and sp_doppler_bin   
    if (X(i)==1 && Y(i)==1)
        a(i)=(A(2,1)+A(1,2))/2;
    elseif (X(i)==1 && Y(i)==7)
        a(i)=(A(2,7)+A(1,6))/2;
    else
        a(i) = interp2(A,Y(i),X(i),'linear');
    end   
end
%disp(['delay = ',num2str(delay(i)),'; doppler(i) = ',num2str(doppler(i)),'; a = ',num2str(a(i))]);

sigma2=(a.*DDMobs(bin_index)').^2+sigma2_n; %1XN
sigma=sqrt(sigma2);
R(logical(eye(N)))=sigma2';

%% Compute off-diagonal values R(N,N)
%sigma_ij=sigma_i*sigma_j*f
for i=1:N-1
    for j=i+1:N
        index1=12*(Y(i)-1)+X(i); %1-84
        index2=12*(Y(j)-1)+X(j); %1-84
        a1=interp2(A1,index1,index2,'linear');
        b1=interp2(B1,index1,index2,'linear');
        c1=interp2(C1,index1,index2,'linear');
        f=a1+b1/sp_ws+c1/(sp_ws^2); %f=0;
        R(i,j)=sigma(i)*sigma(j)*f; %NXN
        R(j,i)=R(i,j);
    end
end

%% recondition method
if (strcmp(method,'RR')) %ridge regression
    lambda=eig(R);
    delta=(max(lambda)-min(lambda)*k_max)/(k_max-1);
    R1=R+delta*eye(N);
elseif (strcmp(method,'ME')) % Minimum eigenvalue
    [V,D] = eig(R);
    lambda=diag(D);
    D1=zeros(N,N);
    T=max(abs(lambda))/k_max; %threshold eigenvalue
    for i=1:N
        if (abs(lambda(i))>T)
            D1(i,i)=abs(lambda(i));
        else
            D1(i,i)=T;
        end               
    end
    R1=V*D1*V';
end

%% compute inverse
iR=pinv(R1);


%% plot R and iR
%{
iR=pinv(R);
disp(['condition = ',num2str(cond(R))]);
figure('units','normalized','outerposition',[0.2 0.4 0.5 0.4]); %left bottom width height
subplot(1,2,1);imagesc(R);colorbar;
xlabel('effective bins');ylabel('effective bins');title('Covariance matrix');
axis equal;axis tight;
set(gca,'fontsize',16);
subplot(1,2,2);imagesc(iR);colorbar;
xlabel('effective bins');ylabel('effective bins');title('Inverse of covariance matrix');
axis equal;axis tight;
set(gca,'fontsize',16);

%% try with diagonal matrix
R_diag=zeros(N,N);
R_diag(logical(eye(N)))=diag(R);
figure;
subplot(1,2,1);imagesc(R_diag);colorbar;
xlabel('effective bins');ylabel('effective bins');title('Covariance matrix');
axis equal;axis tight;
set(gca,'fontsize',16);
subplot(1,2,2);imagesc(pinv(R_diag));colorbar;
xlabel('effective bins');ylabel('effective bins');title('Inverse of covariance matrix');
axis equal;axis tight;
set(gca,'fontsize',16);
%}

end

