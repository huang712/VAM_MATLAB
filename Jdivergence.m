function [Jdiv,gdiv] = Jdivergence(size,uv_ana,uv_bkg,scale_div)
% Compute the cost function of divergence and its gradient
% gdiv(size*size*2,1)

n=size*size*2;
uv=uv_ana-uv_bkg;
u=reshape(uv(1:n/2),[size size]);
v=reshape(uv(n/2+1:n),[size size]);
div=divergence(u,v);
div=div.^2;
Jdiv=scale_div*sum(sum(div));

% Gradient
du=zeros(size,size);
dv=zeros(size,size);

% Ignore the gradient at the edge
for i=3:size-2
    for j=3:size-2
        du(i,j)=(u(i,j)-u(i,j-2)+v(i+1,j-1)-v(i-1,j-1))/2-(u(i,j+2)-u(i,j)+v(i+1,j+1)-v(i-1,j+1))/2;
        dv(i,j)=(u(i-1,j+1)-u(i-1,j-1)+v(i,j)-v(i-2,j))/2-(u(i+1,j+1)-u(i+1,j-1)+v(i+2,j)-v(i,j))/2;
    end
end

du=reshape(du,[size*size,1]);
dv=reshape(dv,[size*size,1]);
gdiv=[du;dv]*scale_div;
end

