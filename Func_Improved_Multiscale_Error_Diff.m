%This function is mainly to process the image by the algorithm prospered by
%Yuk-Hee Chan.And more about that can be found in the artical of 
%'Feature-preserving multiscale error diffusion for digital halftoning' in
%the Journal of Electronic Imaging 13(3), 639¨C645 (July 2004).

function Bm= Func_Improved_Multiscale_Error_Diff( Input_Image )
%FUNC_IMPROVED_MULTISCALE_ERROR_DIFF Summary of this function goes here
%Input_Image is the image of gray-level intensity.
%Bm is the output image.
%Read the input image.
tic;
%Input=imread('hair2.jpg');
%Input_Image=rgb2gray(Input);
[m,n]=size(Input_Image);
%N=2^log2(max(m,n));
N=256;
r=log(N)/log(4);
Input_Image_Adapter=imresize(Input_Image,[N,N]);
Im=double(Input_Image_Adapter)/255;

%To determine whether we should inverse the input image.
Invers_Tag=sum(sum(Im))/N^2;
if Invers_Tag>0.5;
    Im=1-Im;
end

%To cut the image into r+1 level
Em=Im;
Bm=zeros(N,N);
E{1}=1;
%Bm_Mark is mainly to mark the Bm.If Bm(i,j) is set value, then,Bm_Mark(i,j)=1;
Bm_Mark=zeros(N,N);
iterate=1;
while(E{1}>0.5)
%The parameters of the filter H
 H=[1 2 1; 2 0 2 ; 1 2 1];
I{r+1}=Im;
E{r+1}=Em;
B{r+1}=Bm;
%Cut the image into r+1 layer.And the second level has 4 x 4 matrix;
  for i=r+1:-1:2;
     N1=4^(i-1);
     I11=double(I{i}(1:4:N1,1:4:N1));   I12=double(I{i}(1:4:N1,2:4:N1));  I13=double(I{i}(1:4:N1,3:4:N1));   I14=double(I{i}(1:4:N1,4:4:N1));
     E11=double(E{i}(1:4:N1,1:4:N1));   E12=double(E{i}(1:4:N1,2:4:N1));  E13=double(E{i}(1:4:N1,3:4:N1));   E14=double(E{i}(1:4:N1,4:4:N1));
     B11=double(B{i}(1:4:N1,1:4:N1));   B12=double(B{i}(1:4:N1,2:4:N1));  B13=double(B{i}(1:4:N1,3:4:N1));   B14=double(I{i}(1:4:N1,4:4:N1));
     
     I21=double(I{i}(2:4:N1,1:4:N1));   I22=double(I{i}(2:4:N1,2:4:N1));  I23=double(I{i}(2:4:N1,3:4:N1));   I24=double(I{i}(2:4:N1,4:4:N1));
     E21=double(E{i}(2:4:N1,1:4:N1));   E22=double(E{i}(2:4:N1,2:4:N1));  E23=double(E{i}(2:4:N1,3:4:N1));   E24=double(E{i}(2:4:N1,4:4:N1));
     B21=double(B{i}(2:4:N1,1:4:N1));   B22=double(B{i}(2:4:N1,2:4:N1));  B23=double(B{i}(2:4:N1,3:4:N1));   B24=double(I{i}(2:4:N1,4:4:N1));
     
     I31=double(I{i}(3:4:N1,1:4:N1));   I32=double(I{i}(3:4:N1,2:4:N1));  I33=double(I{i}(3:4:N1,3:4:N1));   I34=double(I{i}(3:4:N1,4:4:N1));
     E31=double(E{i}(3:4:N1,1:4:N1));   E32=double(E{i}(3:4:N1,2:4:N1));  E33=double(E{i}(3:4:N1,3:4:N1));   E34=double(E{i}(3:4:N1,4:4:N1));
     B31=double(B{i}(3:4:N1,1:4:N1));   B32=double(B{i}(3:4:N1,2:4:N1));  B33=double(B{i}(3:4:N1,3:4:N1));   B34=double(I{i}(3:4:N1,4:4:N1));
     
     I41=double(I{i}(4:4:N1,1:4:N1));   I42=double(I{i}(4:4:N1,2:4:N1));  I43=double(I{i}(4:4:N1,3:4:N1));   I44=double(I{i}(4:4:N1,4:4:N1));
     E41=double(E{i}(4:4:N1,1:4:N1));   E42=double(E{i}(4:4:N1,2:4:N1));  E43=double(E{i}(4:4:N1,3:4:N1));   E44=double(E{i}(4:4:N1,4:4:N1));
     B41=double(B{i}(4:4:N1,1:4:N1));   B42=double(B{i}(4:4:N1,2:4:N1));  B43=double(B{i}(4:4:N1,3:4:N1));   B44=double(I{i}(4:4:N1,4:4:N1));
     I{i-1}=(I11+I12+I13+I14)+(I21+I22+I23+I24)+(I31+I32+I33+I34)+(I41+I42+I43+I44);
     E{i-1}=(E11+E12+E13+E14)+(E21+E22+E23+E24)+(E31+E32+E33+E34)+(E41+E42+E43+E44);
     B{i-1}=(B11+B12+B13+B14)+(B21+B22+B23+B24)+(B31+B32+B33+B34)+(B41+B42+B43+B44); 
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


%Assume that a white dot should be introduced in the region of interest.
%To choose whether the black dot or the white dot should be introduced in the region of interest.
White_Black=1;
%If average of the Em elements in the region is larger than 0.5, it
%indicates that this is a bright region and a black dot should be
%introduced.
if sum(sum(Im(s:s+3,t:t+3)))/16>0.5
    White_Black=0;
end
%To add a 1-pixe frame of value 0 and a 1-pixel frame of value 1 to E and B
%first.
Em_Frame=zeros(N+2,N+2);
Em_Frame(2:N+1,2:N+1)=Em;

Bm_Frame=ones(N+2,N+2);
Bm_Frame(2:N+1,2:N+1)=Bm;

Bm_Mark_Frame=zeros(N+2,N+2);
Bm_Mark_Frame(2:N+1,2:N+1)=Bm_Mark;

%The default type of the dot should be introduced is the same with the
%actual circumstance.
if White_Black==1;
%Updata the error Image E by the method proposed by Chun
Bm_Frame(s_exact+1,t_exact+1)=1;
e=1-Em_Frame(s_exact+1,t_exact+1);
su=sum(sum(H));
%%Adjust the parmaters of the filter of H
for p=1:3;
    for q=1:3;
        if Bm_Mark_Frame(s_exact+p-1,t_exact+q-1)==1;
             H(p,q)=0;
             su=su-1;
        end
    end
end
su=sum(sum(H));
H1=H/su;
for mm=1:3;
    for nn=1:3;
       if Bm_Mark_Frame(s_exact+mm-1,t_exact+nn-1)~=1;
       Em_Frame(s_exact+mm-1,t_exact+nn-1)=Em_Frame(s_exact+mm-1,t_exact+nn-1)-e*H1(mm,nn);
       end
    end
end
%To mark the dot has been assigned.
Em_Frame(s_exact+1,t_exact+1)=0;
Bm_Mark_Frame(s_exact+1,t_exact+1)=1;
%The black dot should introduced in actual.
else 
 
%To compute the total black dots of the region;
Thresh=16-round(sum(sum(Im(s:s+3,t:t+3))));
T=16-sum(sum(Bm_Mark_Frame(s+1:s+4,t+1:t+4)));
if T>Thresh; 
    Bm_Mark_Frame(s+1:s+4,t+1:t+4)=1;  
    Bm_Frame(s+1:s+4,t+1:t+4)=1;  
end
  %Bm_Frame(s+1:s+4,t+1:t+4)=1; 
%Inverse the elements.
   for i=1:4;
        for j=1:4;
            if Bm_Mark_Frame(s+i,t+j)~=1;
            Em_Frame(s+i,t+j)=1-Em_Frame(s+i,t+j);
            Bm_Frame(s+1:s+4,t+1:t+4)=1;
            end
       end
    end
Em=Em_Frame(2:N+1,2:N+1);
 %%Find the maximum error intensity.
for ii=1:3
    for jj=1:3
        Overlapping_Mat(ii,jj)=sum(sum(Em(s+ii-1:s+ii,t+jj-1:t+jj)));
    end
end
[uu,vv]=find(Overlapping_Mat==max(max(Overlapping_Mat)));
MM_UV=max(size(uu));
Rand_Number_U=MM_UV*rand();
Rand_Number_V=MM_UV*rand();
uu=uu(ceil(Rand_Number_U),1);
vv=vv(ceil(Rand_Number_V),1);
%Secondly, we locate the exact pixel which needs to be processed
[uuu,vvv]=find(Em(s+uu-1:s+uu,t+vv-1:t+vv)==max(max(Em(s+uu-1:s+uu,t+vv-1:t+vv))));
MM_UV=max(size(uuu));
Rand_Number_U=MM_UV*rand();
Rand_Number_V=MM_UV*rand();
uuu=uuu(ceil(Rand_Number_U),1);
vvv=vvv(ceil(Rand_Number_V),1);
s_exact=s+uu-1+rem(uuu+1,2);
t_exact=t+vv-1+rem(vvv+1,2);

%initial the first value 
e=1-Em(s_exact,t_exact);
Bm(s_exact,t_exact)=1;
su=sum(sum(H));
for p=1:3
    for q=1:3
        if Bm_Mark_Frame(s_exact+p-1,t_exact+q-1)==1;
             H(p,q)=0;
             su=su-1;
        end
    end
end
su=sum(sum(H));
H1=H/su;
Bm_Mark_Frame(s_exact+1,t_exact+1)=1;
%Error diffusion 
for mm=1:3;
    for nn=1:3;
       if Bm_Mark_Frame(s_exact+mm-1,t_exact+nn-1)~=1;
       Em_Frame(s_exact+mm-1,t_exact+nn-1)=Em_Frame(s_exact+mm-1,t_exact+nn-1)-e*H1(mm,nn);
       end
    end
end
Em_Frame(s_exact+1,t_exact+1)=0;
Bm_Frame(s_exact+1,t_exact+1)=0;
Bm_Mark_Frame(s_exact+1,t_exact+1)=1;

%Inverse the matrix Em
for mm=1:4;
    for nn=1:4;
       if Bm_Mark_Frame(s+mm,t+nn)~=1;
       Em_Frame(s+mm,t+nn)=1-Em_Frame(s+mm,t+nn);
       end
    end
end
end
Em=Em_Frame(2:N+1,2:N+1);
Bm=Bm_Frame(2:N+1,2:N+1);
%fprintf('The s_exact in the iterate is  %d.\n',s_exact);
%fprintf('The t_exact in the iterate is  %d.\n',t_exact);
fprintf('The iterate is  %d,please wait for a moment.\n',iterate);
iterate=iterate+1;
%if iterate>10000;
 % break;
%end
end

%Inverse the ouput again
if Invers_Tag>0.5;
    Bm=1-Bm;
end
toc;



























