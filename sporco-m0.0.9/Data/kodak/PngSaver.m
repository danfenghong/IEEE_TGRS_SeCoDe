in_name = 'kiel.grey.bmp';    %the name of your input image

out_name = 'kiel.grey.png';   %the name of the desired output

IM = imread(in_name);   %read in the image

imwrite(IM, out_name);  %write it out