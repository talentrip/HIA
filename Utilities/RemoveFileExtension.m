function filename = RemoveFileExtension(filename)

% make sure no filename extensions are included in filename
filename(find(filename=='.'):end) = [];