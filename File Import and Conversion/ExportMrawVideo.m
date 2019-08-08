function ExportMrawVideo(filename,ImageMatrix,ImageInfo)

% Remove .mat or .mraw or any file extension
filename = RemoveFileExtension(filename);

% this writes out the camera data in bitN form and spits it out as 16-bits
testmraw_filename = [filename '.mraw'];
fid = fopen(testmraw_filename,'w','b');    % file is saved as a big-endian, hence the 'b'

% One loop iteration per video frame
convertstring = ['ubit' num2str(ImageInfo.BitDepth)];
P = size(ImageMatrix,3);
for i=1:P
    if mod(i,round(P/10)) == 0; disp(num2str(i)); end
    x = ImageMatrix(:,:,i)';
    fwrite(fid, x, convertstring);
end
fclose(fid);