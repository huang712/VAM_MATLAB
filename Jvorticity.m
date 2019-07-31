function [Jvor,gvor] = Jvorticity(size,uv_ana,uv_bkg,scale_vor)
%compute the cost function of vorticity and its gradient
%gvor(size*size*2,1)

n=size*size*2;
uv=uv_ana-uv_bkg;
u=reshape(uv(1:n/2),[size size]);
v=reshape(uv(n/2+1:n),[size size]);
vor=curl(u,v);
vor=vor.^2;
Jvor=scale_vor*sum(sum(vor));

%gradient
du=zeros(size,size);
dv=zeros(size,size);

for i=3:size-2
    for j=3:size-2
        du(i,j)=-(v(i-1,j+1)-v(i-1,j-1)-u(i,j)+u(i-2,j))/8+(v(i+1,j+1)-v(i+1,j-1)-u(i+2,j)+u(i,j))/8;
        dv(i,j)=(v(i,j)-v(i,j-2)-u(i+1,j-1)+u(i-1,j-1))/8-(v(i,j+2)-v(i,j)-u(i+1,j+1)+u(i-1,j+1))/8;
    end
end

du=reshape(du,[size*size,1]);
dv=reshape(dv,[size*size,1]);
gvor=[du;dv]*scale_vor;
end

