%This is a test example for the function of Func_Multiscale_Error_Diff(Input_Image);
Input=imread('hair2.jpg');
Input=rgb2gray(Input);
%You had better resize the image into 256 x 256��or it cost much time to get
%the result.
%Bm is the image after processsing by the algorithm prospered by Chun.
N=256;
Input=imresize(Input,[N,N]);
Bm=Func_Multilevel_Halftoing(Input);
figure(1);
imshow(Bm);
title('The image after processing by the Multilevel Halftoing');