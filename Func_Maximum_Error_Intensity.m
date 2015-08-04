function [s,t]=Func_Maximum_Error_Intensity(Input_Image)
E{1}=1.0;
Em=Input_Image;
[m,n]=size(Input_Image);
N=2^(log2(max(m,n)));
r=log(N)/log(4);
E{r+1}=Em;
%Cut the image into r+1 layer.And the second level has 4 x 4 matrix;
  for i=r+1:-1:2;
     N1=4^(i-1);
     E11=double(E{i}(1:4:N1,1:4:N1));
     E12=double(E{i}(1:4:N1,2:4:N1)); 
     E13=double(E{i}(1:4:N1,3:4:N1)); 
     E14=double(E{i}(1:4:N1,4:4:N1));
   
     E21=double(E{i}(2:4:N1,1:4:N1));  
     E22=double(E{i}(2:4:N1,2:4:N1)); 
     E23=double(E{i}(2:4:N1,3:4:N1));  
     E24=double(E{i}(2:4:N1,4:4:N1));
    
     E31=double(E{i}(3:4:N1,1:4:N1)); 
     E32=double(E{i}(3:4:N1,2:4:N1));
     E33=double(E{i}(3:4:N1,3:4:N1));
     E34=double(E{i}(3:4:N1,4:4:N1));
     
     E41=double(E{i}(4:4:N1,1:4:N1));
     E42=double(E{i}(4:4:N1,2:4:N1));
     E43=double(E{i}(4:4:N1,3:4:N1)); 
     E44=double(E{i}(4:4:N1,4:4:N1));
     E{i-1}=(E11+E12+E13+E14)+(E21+E22+E23+E24)+(E31+E32+E33+E34)+(E41+E42+E43+E44);
  end
  
%Find the first pixel in the 4 x 4 matrix using the maximum error intensity guidance
u=1;
v=1;
i=2;
s=1;
t=1;
% The finding process starts at level 2th
%(s,t) is the absolute coordinates of the first row and first column of the
%4 x 4 matrix.
while(i<=r+1)
    s_1=s;
    t_1=t;
%Get the exact loaction of the loaction firstly.
%We use the overlapping subregions to locate it.
%First,we cut the 4 x 4 matrixe into nine 2 x 2 matrix which is overlapped.
%Find the maximum error intensity. 
 for ii=1:3;
    for jj=1:3;
        Overlapping_Mat(ii,jj)=sum(sum(E{i}(s+ii-1:s+ii,t+jj-1:t+jj)));
    end
end
[uu,vv]=find(Overlapping_Mat==max(max(Overlapping_Mat)));
MM_UV=max(size(uu));
Rand_Number_U=MM_UV*rand();
Rand_Number_V=MM_UV*rand();
uu=uu(ceil(Rand_Number_U),1);
vv=vv(ceil(Rand_Number_V),1);
%Secondly, we locate the exact pixel which needs to be processed
[uuu,vvv]=find(E{i}(s+uu-1:s+uu,t+vv-1:t+vv)==max(max(E{i}(s+uu-1:s+uu,t+vv-1:t+vv))));
MM_UV=max(size(uuu));
Rand_Number_U=MM_UV*rand();
Rand_Number_V=MM_UV*rand();
uuu=uuu(ceil(Rand_Number_U),1);
vvv=vvv(ceil(Rand_Number_V),1);
%(s_exact,t_exact) is the absolute coordinates of the pixel to be
%processed.
s_exact=s+uu-1+rem(uuu+1,2);
t_exact=t+vv-1+rem(vvv+1,2);
    if i==r+1;
        break;
    end
    i=i+1;
    s=4*s_exact-3;
    t=4*t_exact-3;
end 