function  Maxpixelposition = ComputeSpokeposition(spokeSurface,numcumpute)
% ���ú���ComputeSpokeposition������Ƭ�������������ֵ��λ�ã���λΪ�Ƕȣ���Χ
% ��-180��180
% numcumputeΪ��Ҫ�������Ƭ������Ĭ�ϴӵ�һ�ſ�ʼ
Maxpixelposition = zeros(1,numcumpute);
for i=1:numcumpute
   [xposition yposition] = find(spokeSurface(i,:)==max(spokeSurface(i,:)));
   Maxpixelposition(i) = yposition*2-180;
end
