function  Maxpixelposition = ComputeSpokeposition(spokeSurface,numcumpute)
% 调用函数ComputeSpokeposition计算照片矩阵中最大像素值的位置，单位为角度，范围
% 从-180到180
% numcumpute为想要计算的照片张数，默认从第一张开始
Maxpixelposition = zeros(1,numcumpute);
for i=1:numcumpute
   [xposition yposition] = find(spokeSurface(i,:)==max(spokeSurface(i,:)));
   Maxpixelposition(i) = yposition*2-180;
end
