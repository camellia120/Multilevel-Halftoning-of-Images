% we use the quadtree method to find the best needed pixel to scan
% And this function is mainly to process the image by the multiscal error
% diffusion prospered by the Giorgos Sarailidis and Ioannis Katsavounidis.
% To learn more about the algorithm, you can read the artical of 'A Multiscale Error Diffusion
%Technique for Digital Halftoning' by Ioannis Katsavounidis and C.-C. Jay
%Kuo.And u can find the artical in IEEE TRANSACTIONS ON IMAGE PROCESSING, VOL. 6, NO. 3, MARCH 1997

%Input:Input_Image,which should be fristly read by the function of imread;
%Output:Output_Image is the image after processing by the algorithm
%And the pixel of the output is of value 0 or 1;


function Output_Image=Func_Multiscale_Error_Diff(Input_Image)
tic;
%when the size of input image is 256 x 256,the iteration is 30195 and 
%the elapsed time is 77.681618 seconds.And the size if the input image 
%is 128 x 128,the elapsed time 39.9657 seconds,so we set the size of the
%image to 256 x 256
%to expand the image to own the size of N x N,where N=2^r
%Input_Image1=imread('hair2.jpg');
%Input_Image=rgb2gray(Input_Image1);
[m,n] = size(Input_Image);
N=2^(floor(log2(max(m,n))));
Input_Image_Adapter=imresize(Input_Image,[N,N]);
%to show the image after resizing
%figure(1);
%imshow(Input_Image_Adapter);
r=log2(N);
Im=double(Input_Image_Adapter)/255;
%The output image is regarded as four scal gray. B is limited to the set of {0,1/3,2/3,1}
%Initial the condition
%You can choose another method to initial the Em and Bm
Bm=zeros(N,N);
Em=Im;
%And Im,Em,Bm stands for the level (r+1)th,respectively,it equals to the 
%use the method of maximum intensity method to find the best needed pixel
%cut the input image into r+1 layers
%initial the first value of error and output image
%iterate denotes the iteration
iterate=1;
%when the E{1}is less than 0.5,the algorithm stops.
%The value of the last level is the value of the input image 
E{1}=1.0;
while(E{1}>=0.5)
I{r+1}=Im;
B{r+1}=Bm;
E{r+1}=Em;
for i=r+1:-1:2;
     N1=2^(i-1);
     I1=double(I{i}(1:2:N1,1:2:N1));
     E1=double(E{i}(1:2:N1,1:2:N1));
     B1=double(B{i}(1:2:N1,1:2:N1));
     
     I2=double(I{i}(1:2:N1,2:2:N1));
     E2=double(E{i}(1:2:N1,2:2:N1));
     B2=double(B{i}(1:2:N1,2:2:N1));
     
     I3=double(I{i}(2:2:N1,1:2:N1));
     E3=double(E{i}(2:2:N1,1:2:N1));
     B3=double(B{i}(2:2:N1,1:2:N1));
     
     I4=double(I{i}(2:2:N1,2:2:N1));
     E4=double(E{i}(2:2:N1,2:2:N1));
     B4=double(B{i}(2:2:N1,2:2:N1));
     
     I{i-1}=(I1+I2+I3+I4);
     E{i-1}=(E1+E2+E3+E4);
     B{i-1}=(B1+B2+B3+B4);    
end

%After that, the input image I,the output image B,the error E have been
%divided into many layers

%Find the needed pixel as the maximum intensity guidance
u=1;
v=1;
i=2;
s=1;
t=1;
% The finding process starts at level 2th
while(i<=r+1)
    s_1=s;
    t_1=t;
    [u,v]=find(E{i}(s:s+1,t:t+1)==max(max(E{i}(s:s+1,t:t+1))));
    MM=max(size(u));
    rand_number1=MM*rand();
    rand_number2=MM*rand();
    u=u(ceil(rand_number1),1);
    v=v(ceil(rand_number2),1);
    if i==r+1;
        break;
    end
    i=i+1; 
    if s_1*t_1*u*v==1
        s=1;
        t=1;
        continue;
    end
    if u==1;
        s=2*s-1;
    else if u==2;
         s=2*s+1;
        end
    end
    if v==1;
        t=2*t-1;
    else if v==2;
        t=2*t+1;
        end
    end
end   
% the best maximum error intensity pixel
% to revise the error
 s=s+rem(u+1,2);
 t=t+rem(v+1,2);
  %fprintf('this fact is that s=%d.\n',s);
  % fprintf('this fact is that t=%d.\n',t);
e=Em(s,t)-1;
Em(s,t)=e;
Bm(s,t)=1;
%You can also choose this method to revise the E and B,but we recommemd
%the method above
%for k=4:-1:2;
%  if (Threshold(k)<=T)||(Threshold(k-1)<=Bm(s,t));
%      Em(s,t)=T-Threshold(k);
%      Bm(s,t)=Threshold(k);
%      break;
%   end
%end
% use the filter to refresh the error 
% the parameters of noncausal diffusion filter which is refered to the references.

% The parmeters of the inner filter 
if Em(s,t)<0;
Filter_Inner=(1/12)*[1 2 1;2 -12 2;1 2 1];

%The parmeters of the boundary filter
F11=1/5*[-5 2;2 1];
F12=1/5*[2 -5;1 2];
F21=1/5*[2 1;-5 2];
F22=1/5*[2 1;2 -5];
F1_2=1/8*[2 -8 2;1 2 1];
F2_2=1/8*[1 2 1;2 -8 1];
F1_1=1/8*[2 1;-8 2;2 1];
F2_1=1/8*[1 2;2 -8;1 2];
%Process the inner pixel of the image
if s>1&&s<N&&t>1&&t<N-1;
Em(s-1:s+1,t-1:t+1)=Em(s-1:s+1,t-1:t+1)+e*Filter_Inner;

%Process the boundary of the image
else if s*t==1;
       Em(1:2,1:2)=Em(1:2,1:2)+e*F11;
    else if s==1&&t==N;
       Em(1:2,N-1:N)=Em(1:2,N-1:N)+e*F12;
        else if s==N&&t==1;
        Em(N-1:N,1:2)=Em(N-1:N,1:2)+e*F21;
            else if s*t==N^2;
        Em(N-1:N,N-1:N)=Em(N-1:N,N-1:N)+e*F22;
                else if s==1;
        Em(1:2,t-1:t+1)=Em(1:2,t-1:t+1)+e*F1_2;
                    else if s==N;
        Em(N-1:N,t-1:t+1)=Em(N-1:N,t-1:t+1)+e*F2_2;
                        else if t==1;
        Em(s-1:s+1,1:2)=Em(s-1:s+1,1:2)+e*F1_1;
                            else if t==N;
        Em(s-1:s+1,N-1:N)=Em(s-1:s+1,N-1:N)+e*F2_1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
end
fprintf('this will be  %d iterations.\n',iterate);
fprintf('this will be stop when %d is less than 0.5.\n',E{1});
iterate=iterate+1;
%if iterate>500;
%break;
%end
end 
%Adapt the value override.
%save('Bm_Multiscale.mat','Bm');
%save('Em_Multiscale.mat','Em');
Output_Image=Bm;
%figure(2);
%imshow(Bm);
toc;

            
                
            




