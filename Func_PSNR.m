%This function is mainly to evaluate the image after process

%The Input_Image and the Output_Image should be the same size.
%And the color image should firstly be transfered to the gray image.
%Then we calculate the Mean Square Error (MSE) and the Peak Signal to Noise
%Ratio(PSNR)
function PSNR = Func_PSNR(Input_Image,Output_Image)
[m,n]=size(Input_Image);
for i=1:m;
    for j=1:n
        MSE_Element(i,j)=(Input_Image(i,j)-Output_Image(i,j))^2;
    end
end
MSE=double(sum(sum(MSE_Element))/(m*n));
Max_Gray=double(max(max(Input_Image)));
PSNR=10*log10(Max_Gray^2/MSE);