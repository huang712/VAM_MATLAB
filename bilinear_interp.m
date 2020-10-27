function [bi_index, bi_weight] = bilinear_interp(x_vec,y_vec,x,y,resolution)
% Given the lat_vec, lon_vec, and lat/lon of a point, conpute the bi_index[4] and bi_weight[4] of the four points  
% x=lon, y=lat, index=(y-1)*size_x+x
size_x=length(x_vec);
size_y=length(y_vec);

% Find ix1 and ix2, two index that are nearest to x
for i=1:size_x
    if(abs(x_vec(i)-x)<resolution)
        ix1=i;
        ix2=i+1;
        x1=x_vec(ix1);
        x2=x_vec(ix2);
        break;
    end
end

for i=1:size_y
    if(abs(y_vec(i)-y)<resolution)
        iy1=i;
        iy2=i+1;
        y1=y_vec(iy1);
        y2=y_vec(iy2);
        break;
    end
end

bi_index(1) = (iy1-1)*size_x+ix1; % index of Q11
bi_index(2) = (iy1-1)*size_x+ix2; % index of Q21
bi_index(3) = (iy2-1)*size_x+ix1; % index of Q12
bi_index(4) = (iy2-1)*size_x+ix2; % index of Q22

bi_weight(1) = (x2-x)*(y2-y)/((x2-x1)*(y2-y1)); % weight of Q11
bi_weight(2) = (x-x1)*(y2-y)/((x2-x1)*(y2-y1)); % weight of Q21
bi_weight(3) = (x2-x)*(y-y1)/((x2-x1)*(y2-y1)); % weight of Q12
bi_weight(4) = (x-x1)*(y-y1)/((x2-x1)*(y2-y1)); % weight of Q22

end

