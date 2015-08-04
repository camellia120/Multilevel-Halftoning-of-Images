%This function is mainly to process the image by the algorithm prospered by
%Yuk-Hee Chan.And more about that can be found in the artical of 
%'Feature-preserving multiscale error diffusion for digital halftoning' in
%the Journal of Electronic Imaging 13(3), 639¨C645 (July 2004).

%FUNC_IMPROVED_MULTISCALE_ERROR_DIFF Summary of this function goes here
%Input_Image is the image of gray-level intensity.
%Bm is the output image.
%Read the input image.
%function Bm=Func_Improved_Multiscale_Error_Diff(Input)
tic;
Input=imread('ramp.png');
Input_Image=rgb2gray(Input);
%Input_Image=Input;
%[m,n]=size(Input_Image);
%N=2^log2(max(m,n));
N=256;
r=log(N)/log(4);
Input_Image_Adapter=imresize(Input_Image,[N,N]);
Im=double(Input_Image_Adapter)/255;
figure(1);
imshow(Im);
%To determine whether we should inverse the input image.
Invers_Tag=sum(sum(Im))/(N^2);
if Invers_Tag<0.5;
    Im=1-Im;
end 

%To cut the image into r+1 level
Im_L=zeros(N+12,N+12);
Im_L(7:N+6,7:N+6)=Im;
figure(2);
imshow(Im_L);
%Em=Im;
Bm_L=zeros(N+12,N+12)+2;
E{1}=1;
%Bm_Mark is mainly to mark the Bm.If Bm(i,j) is set value, then,Bm_Mark(i,j)=1;
Bm_Mark_L=ones(N+12,N+12);
Bm_Mark_L(7:N+6,7:N+6)=0;
Bm_Mark_M{1}=(N+12)^2;
iterate=1;
%To add a 1-pixe frame of value 0 and a 1-pixel frame of value 1 to E and B
%first.
%Em_Frame=zeros(N+2,N+2);
%Em_Frame(2:N+1,2:N+1)=Em;

%Bm_Frame=ones(N+2,N+2);
%Bm_Frame(2:N+1,2:N+1)=Bm;

%Bm_Mark_Frame=ones(N+2,N+2);
%Bm_Mark_Frame(2:N+1,2:N+1)=Bm_Mark; % to edit
%Total number of black dots.
Thresh=N^2-sum(sum(Im));
Black_Dots=0;
while(E{1}>0.5);
 %A simple trick to eliminate the boundary effect.
 
 x0=randi([4,10],1);
 y0=randi([4,10],1);
 Em=Im_L(x0:x0+255,y0:y0+255);
 Bm_Mark=Bm_Mark_L(x0:x0+255,y0:y0+255);
 Bm=Bm_L(x0:x0+255,y0:y0+255);
 figure(3);
 imshow(Em);
 %x_offset=ceil(3*rand());
 %y_offset=ceil(3*rand());
 %Em=Em_Frame(x_offset:x_offset+255,y_offset:y_offset+255);
E{r+1}=Em; 
Bm_Mark_M{r+1}=1-Bm_Mark; 

%Cut the image into r+1 layers
%1 to 4x4 mapping
  for i=r+1:-1:2;
     N1=4^(i-1);
     E11=double(E{i}(1:4:N1,1:4:N1)); 
     E12=double(E{i}(1:4:N1,2:4:N1));     
     E13=double(E{i}(1:4:N1,3:4:N1)); 
     E14=double(E{i}(1:4:N1,4:4:N1));
     
     Bm_Mark11=Bm_Mark_M{i}(1:4:N1,1:4:N1);
     Bm_Mark12=Bm_Mark_M{i}(1:4:N1,2:4:N1);
     Bm_Mark13=Bm_Mark_M{i}(1:4:N1,3:4:N1);
     Bm_Mark14=Bm_Mark_M{i}(1:4:N1,4:4:N1);
      
     E21=double(E{i}(2:4:N1,1:4:N1)); 
     E22=double(E{i}(2:4:N1,2:4:N1));  
     E23=double(E{i}(2:4:N1,3:4:N1));  
     E24=double(E{i}(2:4:N1,4:4:N1));
     
     Bm_Mark21=Bm_Mark_M{i}(2:4:N1,1:4:N1);
     Bm_Mark22=Bm_Mark_M{i}(2:4:N1,2:4:N1);
     Bm_Mark23=Bm_Mark_M{i}(2:4:N1,3:4:N1);
     Bm_Mark24=Bm_Mark_M{i}(2:4:N1,4:4:N1);

     E31=double(E{i}(3:4:N1,1:4:N1)); 
     E32=double(E{i}(3:4:N1,2:4:N1)); 
     E33=double(E{i}(3:4:N1,3:4:N1));  
     E34=double(E{i}(3:4:N1,4:4:N1));
     
     Bm_Mark31=Bm_Mark_M{i}(3:4:N1,1:4:N1);
     Bm_Mark32=Bm_Mark_M{i}(3:4:N1,2:4:N1);
     Bm_Mark33=Bm_Mark_M{i}(3:4:N1,3:4:N1);
     Bm_Mark34=Bm_Mark_M{i}(3:4:N1,4:4:N1);
    
     E41=double(E{i}(4:4:N1,1:4:N1));  
     E42=double(E{i}(4:4:N1,2:4:N1));  
     E43=double(E{i}(4:4:N1,3:4:N1)); 
     E44=double(E{i}(4:4:N1,4:4:N1));
     
     Bm_Mark41=Bm_Mark_M{i}(4:4:N1,1:4:N1);
     Bm_Mark42=Bm_Mark_M{i}(4:4:N1,2:4:N1);
     Bm_Mark43=Bm_Mark_M{i}(4:4:N1,3:4:N1);
     Bm_Mark44=Bm_Mark_M{i}(4:4:N1,4:4:N1);
     
     E{i-1}=(E11+E12+E13+E14)+(E21+E22+E23+E24)+(E31+E32+E33+E34)+(E41+E42+E43+E44);
     Bm_Mark_M{i-1}=(Bm_Mark11+Bm_Mark12+Bm_Mark13+Bm_Mark14)+(Bm_Mark21+Bm_Mark22+Bm_Mark23+Bm_Mark24)+(Bm_Mark31+Bm_Mark32+Bm_Mark33+Bm_Mark34)+(Bm_Mark41+Bm_Mark42+Bm_Mark43+Bm_Mark44);
  end
  
%Find the first pixel in the 4 x 4 matrix using the maximum error intensity guidance
u=1; %local coordinates in 3x3
v=1;
i=2; %start with layer no.2, which is 4x4
s=1; %exact coordinates in 4x4
t=1;

%To compute the total black dots of the region;
while(i<=r+1)
    Inverse=0;
    s_1=s;
    t_1=t;
%Get the u,v loaction firstly.
%We use the overlapping subregions to locate it.
%First,we cut the 4 x 4 matrixe into nine 2 x 2 matrix which is overlapped.
%Find the maximum error intensity one out of nine 2 x 2. 
Max=-100000;
uu_vv=[];
 for ii=1:3;
    for jj=1:3;
        Overlapping_Mat(ii,jj)=sum(sum(E{i}(s+ii-1:s+ii,t+jj-1:t+jj)));
        Total_Number(ii,jj)=sum(sum(Bm_Mark_M{i}(s+ii-1:s+ii,t+jj-1:t+jj)));
        if Total_Number(ii,jj)==0; % all are assigned and marked
          
            continue;
        else
            if Overlapping_Mat(ii,jj)>=Max;
                Max=Overlapping_Mat(ii,jj);
                uu_vv=[uu_vv;ii,jj];
            end            
        end
    end
 end
uu=uu_vv(end,1);
vv=uu_vv(end,2);
%Secondly, we locate the exact element one out of 2 x 2
uuu_vvv=[];
Max=-100000;
for ii=1:2;
    for jj=1:2;
        Total_Number(ii,jj)=Bm_Mark_M{i}(s+uu+ii-2,t+vv+jj-2);
        if Total_Number(ii,jj)==0; % all are assigned and marked
            continue;
        else
            if E{i}(s+uu+ii-2,t+vv+jj-2)>=Max;
                Max=E{i}(s+uu+ii-2,t+vv+jj-2);
                uuu_vvv=[uuu_vvv;ii,jj];
            end 
        end
    end
end
uuu=uuu_vvv(end,1);
vvv=uuu_vvv(end,2);
%(s_exact,t_exact) is the absolute coordinates of the pixel to be
%processed.
s_exact=s+uu-1+rem(uuu+1,2);
t_exact=t+vv-1+rem(vvv+1,2); 
i=i+1;
     if i==r+2;
        break;
     end
     s=4*s_exact-3;
     t=4*t_exact-3;  
     
%%%To decide whether we shoulde inverse the matrix.
     if i==4;% 4 or 5 test
         if E{i-1}(s_exact,t_exact)/Bm_Mark_M{i-1}(s_exact,t_exact)>0.5; %find a black dot in bright region 
             Inverse=1;
         end
     end
     
if Inverse==1;
     while(i<=r+1);
  
%Find the minimum error intensity one out of nine 2 x 2. 
Min=100000;
uu_vv=[];
 for ii=1:3;
    for jj=1:3;
        Overlapping_Mat(ii,jj)=sum(sum(E{i}(s+ii-1:s+ii,t+jj-1:t+jj)));
        Total_Number(ii,jj)=sum(sum(Bm_Mark_M{i}(s+ii-1:s+ii,t+jj-1:t+jj)));
        if Total_Number(ii,jj)==0;
            continue;
        else
            if Overlapping_Mat(ii,jj)<=Min;
                Min=Overlapping_Mat(ii,jj);
                uu_vv=[uu_vv;ii,jj];
            end
            
        end
    end
 end
uu=uu_vv(end,1);
vv=uu_vv(end,2);
%Secondly, we locate the exact element ont out of 2 x 2
uuu_vvv=[];
Min=100000;
for ii=1:2;
    for jj=1:2;
        Total_Number(ii,jj)=Bm_Mark_M{i}(s+uu+ii-2,t+vv+jj-2);
        if Total_Number(ii,jj)==0;
            continue;
        else
            if E{i}(s+uu+ii-2,t+vv+jj-2)<=Min;
                Min=E{i}(s+uu+ii-2,t+vv+jj-2);
                uuu_vvv=[uuu_vvv;ii,jj];
            end 
        end
    end
end
uuu=uuu_vvv(end,1);
vvv=uuu_vvv(end,2);
%(s_exact,t_exact) is the absolute coordinates of the pixel to be
%processed.
s_exact=s+uu-1+rem(uuu+1,2);
t_exact=t+vv-1+rem(vvv+1,2); 
i=i+1;
     if i==r+2;
        break;
     end
     s=4*s_exact-3;
     t=4*t_exact-3; 
   end
end

 %Refresh the layers.
 end 
%Assume that by default a white dot should be introduced in the region of interest.
%Actually whether black dot or white dot should be assigned depends on whehter the region has been inverted 
if Inverse==0;
Bm_L(x0+s_exact-1,y0+t_exact-1)=1;
%Updata the error Image E by the method proposed by Chun
e=1-Em(s_exact,t_exact);
%To set the corresponding dot in Em_Frame to zero.
Im_L(x0+s_exact-1,y0+t_exact-1)=0;
else
Bm_L(x0+s_exact-1,y0+t_exact-1)=0;
e=-Em(s_exact,t_exact);
Im_L(x0+s_exact-1,y0+t_exact-1)=0;
Black_Dots=Black_Dots+1;
end

Bm_Mark_L(x0+s_exact-1,y0+t_exact-1)=1;

%The parameters of the filter H

d=1;
su=0;
%%Adjust the parmaters of the filter of H
while(su==0);
     for p=1:2*d+1;
         for q=1:2*d+1;
             
                 if Bm_Mark_L(x0+s_exact-1+p-d,y0+t_exact-1+q-d)==1;
                    H(p,q)=0;
                 else
                    H(p,q)=2*d+1-(abs(p-d-1)+abs(q-d-1));  
                 end
            
         end
     end
    su=sum(sum(H));
    if su==0;
       d=d+1;
    end
    
     if d==3;         
       break;
     end
end

if su==0;
   e_equal=e/Bm_Mark_M{1};
   Em=Em-Bm_Mark_M{r+1}*e_equal;
   Im_L(x0:x0+255,y0:y0+255)=Em;
else

% Change the filter H 
H1=H/su;
%To diffuse the error with the filter of 3 x 3 matrix H
for mm=1:2*d+1;
    for nn=1:2*d+1;
%We just process the pixels in the Em matrix.
        %if (s_exact+mm-d>=1)&&(s_exact+mm-d<=N)&&(t_exact+nn-d>=1)&&(t_exact+nn-d<=N);
              Im_L(x0+s_exact-1+mm-d,y0+t_exact-1+nn-d)=Im_L(x0+s_exact-1+mm-d,y0+t_exact-1+nn-d)-e*H1(mm,nn);
       % end
    end
end  

end


fprintf('The s_exact in the iterate is  %d.\n',s_exact+x0);
fprintf('The t_exact in the iterate is  %d.\n',t_exact+y0);
fprintf('The iterate is  %d, please wait for a moment.\n',iterate);
iterate=iterate+1;


if Black_Dots>=Thresh;
    break;
end
%if iterate>1000;
%  break;
%end
end

%Inverse the ouput again
if Invers_Tag<0.5; % to edit
    Bm=1-Bm;
end
figure(4);
imshow(Bm);
toc;