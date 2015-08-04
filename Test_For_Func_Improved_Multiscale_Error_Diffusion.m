%This is a test example for the function of Func_Improved_Multiscale_Error_Diff(Input);

Input=imread('hair1.jpg');
Input=rgb2gray(Input);
%You had better resize the image into 256 x 256£¬or it cost much time to get
%the result.
%Bm is the image after processsing by the algorithm prospered by Chun.
N=256;
Input=imresize(Input,[N,N]);
%Bm is the image after processsing by the algorithm prospered by Chun.
Bm=Func_Improved_Multiscale_Error_Diff(Input);
figure(1);
imshow(Bm);
title('The image after processing by the Improved Multiscale Error Diffusion');


