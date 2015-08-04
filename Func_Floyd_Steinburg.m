%Floyd-Steinberg error diffusion.
function Output_Image=Func_Floyd_Steinburg(Input_Image)
img=Input_Image;
[m,n]=size(img);
% re denotes the output image
re=zeros(m,n);
tmp=zeros(m+2,n+2);
tmp(2:m+1,2:n+1)=img;

for i=2:m+1
    for j=2:n+1
        if tmp(i,j)<128
            re(i-1,j-1)=0;
            err=tmp(i,j);
        else
            re(i-1,j-1)=255;
            err=tmp(i,j)-255;
        end
        tmp(i,j+1)=tmp(i,j+1)+7/16*err;
        tmp(i+1,j-1)=tmp(i+1,j-1)+3/16*err;
        tmp(i+1,j)=tmp(i+1,j)+5/16*err;
        tmp(i+1,j+1)=tmp(i+1,j+1)+1/16*err;
    end
end
Output_Image=re;
end