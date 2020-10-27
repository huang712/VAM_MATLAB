function [Jlap,glap] = Jlaplacian(size,uv_ana,uv_bkg,scale_lap)
% Compute the cost function of Laplacian and its gradient

n=size*size*2;
uv=uv_ana-uv_bkg;
u=reshape(uv(1:n/2),[size size]);
v=reshape(uv(n/2+1:n),[size size]);

[ulap,dulap]=laplacian(u);
[vlap,dvlap]=laplacian(v);

Jlap=scale_lap*(ulap+vlap);
du=reshape(dulap,[size*size,1]);
dv=reshape(dvlap,[size*size,1]);
glap=[du;dv]*scale_lap;

end

