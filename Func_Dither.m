%dither
function Output_Image=Func_dither(Input_Image)
dif=dither(Input_Image);
Output_Image=dif;
Output_Image=255*uint8(Output_Image);
end