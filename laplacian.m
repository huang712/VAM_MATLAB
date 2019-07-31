function [Jlap, dlap]= laplacian(A)
    %A is a matrix
    %Use five point stencil method: https://en.wikipedia.org/wiki/Discrete_Laplace_operator
    
    [x,y]=size(A);
    lap1 = zeros(x+2,y+2); %expand lap1
    A1 = zeros(x+2,y+2); %expand A to compute finite difference
    A1(2:x+1,2:y+1)=A;
    A1(1,2:y+1)=A(1,:);
    A1(x+2,2:y+1)=A(x,:);
    A1(2:x+1,1)=A(:,1);
    A1(2:x+1,y+2)=A(:,y);
    
    for i=2:x+1
        for j=2:y+1
            lap1(i,j)=A1(i-1,j)+A1(i+1,j)+A1(i,j-1)+A1(i,j+1)-4*A1(i,j);
        end
    end
    lap=lap1(2:x+1,2:y+1);
    lap=lap.^2;
    Jlap=sum(sum(lap));
    
    dlap=zeros(x,y);
    for i=3:x-2
        for j=3:y-2
            dlap(i,j)=2*(A(i,j-2)+A(i-2,j)+A(i+2,j)+A(i,j+2)...
                +2*A(i-1,j-1)+2*A(i-1,j+1)+2*A(i+1,j-1)+2*A(i+1,j+1)...
                -8*A(i,j-1)-8*A(i-1,j)-8*A(i+1,j)-8*A(i,j+1)...
                +20*A(i,j));
        end
    end
end

