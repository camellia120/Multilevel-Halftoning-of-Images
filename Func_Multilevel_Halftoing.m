%This function is mainly to process the image by the Alogtithm of
%Multilevel Halftoing.
%Note:This function will call the function of Func_Maximum_Error_Intensity
%and Func_Minimum_Error_Intensity.So please ensure the directory has
%included them firstly. 
function B_Final=Func_Multilevel_Halftoing(Input_Image)
Input=Input_Image;
%set the level of gray of the output image.
m=8;
N=256;
%[k,n]=size(Input);
%N=2^(log2(max(m,n)));
Input=double(Input)/255;
Am=imresize(Input,[N,N]);
A0=ones(N,N);
A{1}=A0-(1-Am).^(m-1);
B{1}=Am;
%The input image is first decomposed into m-1 images
for d=2:m-1
    A{d}=A{d-1}-Am.^(d-1).*(factorial(m-1)/(factorial(d-1)*factorial(m-d))).*(1-Am).^(m-d);
    B{d}=B{d-1};
end

%Initialize mask R(x,y)=1 for all (x,y) to indicate they have not been
%processed.
iterate=1;
R=ones(N,N);
for n=1:(m-1)/2;
    N0=sum(sum(R));
    Budget1=round(sum(sum(A{m-n})));
    Budget0=round(N0-(sum(sum(A{n}))));
    Budget_Ratio=Budget1/Budget0;
    while Budget1+Budget0>0;
        if Budget1>=Budget_Ratio*Budget0;
 %Energy plane E=A{m-n}
         E=A{m-n};
 %Search for a pixel loaction in E via Maximum Error Intensity Guidance
 %We adopte the method prospered in Ref12.And we call the function of
 %Func_Maximun_Error_Intensity to get the coordinate in simplicity.
        [s,t]=Func_Maximum_Error_Intensity(E);
              for i=n:m-n;
                  B{i}(s,t)=1;
 %Diffuse error of A{i}(s,t) to its neighbors in A{i}
 %In simplicty,we add frame to A{i}
                 A_Frame=zeros(N+2,N+2);
                 A_Frame(2:N+1,2:N+1)=A{i};
                 R_Frame=ones(N+2,N+2);
                 R_Frame(2:N+1,2:N+1)=R;
                 H=[1,2,1;2 0 2;1 2 1];
%Adapte the parameters of the filter H1.               
                  for u_x=1:3;
                     for v_y=1:3
                         if R_Frame(s+u_x-1,t+v_y-1)==0;
                            H(u_x,v_y)=0;
                         end
                     end
                 end
                 su=sum(sum(H));
                 if su~=0;
                 H1=H/su;
%Diffuse the error to its neighbors.
                 for p=1:3;
                     for q=1:3;
                         if R_Frame(s+p-1,t+q-1)==0;
                             A{i}(s+p-2,t+q-2)=A{i}(s+p-2,t+q-2)-(B{i}(s,t)-A{i}(s+p-2,t+q-2))*H1(p,q);
                         end
                     end
                 end
                 end
                 A{i}(s,t)=0;
              end
        Budget1=Budget1-1;
        else
 %Energy plane E=A{m-n}
        E=A{n};
 %Search for a pixel loaction in E via Minimum Error Intensity Guidance
 %We adopte the method prospered in Ref12.And we call the function of
 %Func_Minimun_Error_Intensity  to get the coordinate in simplicity.
         [s,t]=Func_Minimum_Error_Intensity(E);
              for i=n:m-n;
                  B{i}(s,t)=0;
 %Diffuse error of A{i}(s,t) to its neighbors in A{i}
 %In simplicty,we add frame to A{i}
                 A_Frame=zeros(N+2,N+2);
                 A_Frame(2:N+1,2:N+1)=A{i};
                 R_Frame=ones(N+2,N+2);
                 R_Frame(2:N+1,2:N+1)=R;
%Adapte the parameters of the filter H1.
                 H=[1,2,1;2 0 2;1 2 1];
                 for u_x=1:3;
                     for v_y=1:3
                         if R_Frame(s+u_x-1,t+v_y-1)==0;
                            H(u_x,v_y)=0;
                         end
                     end
                 end
                 su=sum(sum(H));
  %To adapt the parmater to aviod
                 if su~=0;
                 H1=H/su;
                 for p=1:3;
                     for q=1:3;
                         if R_Frame(s+p-1,t+q-1)==0;
                             A{i}(s+p-2,t+q-2)=A{i}(s+p-2,t+q-2)-(B{i}(s,t)-A{i}(s+p-2,t+q-2))*H1(p,q);
                         end
                     end
                 end
                 end
                  A{i}(s,t)=0;
          end
        Budget0=Budget0-1;
        end
        
 %Update mask R(s,t)=0 to indicate pixel(s,t) was processed.
         R(s,t)=0;
    end
 %Adapte the value in B
    for i=1:N;
        for j=1:N;
              if B{n}(i,j)~=0;
                 B{n}(i,j)=1;
              end
              if B{m-n}(i,j)~=1;
                 B{m-n}(i,j)=0;
               end
        end
    end
   fprintf('the iterate is %d.\n',iterate);
  iterate=1+iterate;
    end
  
B_Sum=zeros(N,N);
%Multilevel halftone is B(x,y)=sum(B{i})/m-1;
for i=1:m-1;
    B_Sum=B_Sum+B{i};
end
B_Final=B_Sum/(m-1);
figure(1);
imshow(B_Final);