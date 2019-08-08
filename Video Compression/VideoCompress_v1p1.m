function [Final8BitSparsity,Data] = VideoCompress_v1p1(Data)

% Optional steps for images
imageopt = 0;
if imageopt
    MeanData = mean(Data,3);
    % normalization?
    
    Data = bsxfun(@minus,Data,MeanData);
end


% Haar transform on time-wise dimension is not as effective as simply
% subtracting out mean image ahead of time, and it is much more
% time-consuming in the compression process.  Don't do it.
% if size(Data,3)~=1
%     disp([sprintf('\r') 'Doing 3D Haar' sprintf('\r')])
%     Data = HaarMultiDim(Data,3,ParallelOK);
% end
% toc

MATLABVersion = version('-release')
if str2num(MATLABVersion(1:4)) > 2009
    ParallelOK = 1;
else
        ParallelOK = 0;
end

if ParallelOK
PoolSize = matlabpool('size');
isOpen = PoolSize > 0;
if ~isOpen
    matlabpool open
end
end

Initial8BitSparsity = sum(abs(Data(:))>2^7)/numel(Data);
disp([num2str(100*Initial8BitSparsity) '% of image pixels require > 8 bits before Haar transforms' sprintf('\r')])

tic
% Haar transform rows
disp(['Doing row Haar...' sprintf('\r')])
Data = HaarMultiDim_v1p1(Data,1,ParallelOK);
toc

Mid8BitSparsity = sum(abs(Data(:))>2^7)/numel(Data);
disp([num2str(100*Mid8BitSparsity) '% of image pixels require > 8 bits after row Haar transform' sprintf('\r') ])

% Columns
disp([sprintf('\r') 'Doing column Haar...' sprintf('\r')])
Data = HaarMultiDim_v1p1(Data,2,ParallelOK);
toc

Final8BitSparsity = sum(abs(Data(:))>2^7)/numel(Data);
disp([num2str(100*Final8BitSparsity) '% of image pixels require > 8 bits after row and column Haar transform' sprintf('\r') ])

function input = HaarMultiDim_v1p1(input,dim,ParallelOK)

% This function will take in a vector and do a Haar transform on the
% desired dimension

%% Check that inputs make sense

S = size(input);
if length(S)==2
    S(3) = 1;
end
% Dimension should be between 1 and 3
if dim<1 || dim>3
    ME = MException('VerifyInput:Dimensionality', ...
        'Selected dimension must be between 1 and 3.');
    throw(ME)
end

% Dimension should be less than of equal to # of dimensions of matrix
if dim>length(S)
    ME = MException('VerifyInput:Dimensionality', ...
        'Selected dimension greater than number of matrix dimensions.');
    throw(ME)
end

%% Do the Haar transform on this dimension

if ParallelOK
PoolSize = matlabpool('size');
else
    PoolSize = 1;
end
if dim==1
    
    N = numel(input(1,:,1));
    pad = mod(N,2);
    Q = nextpow2(N);
    X = round([N 2.^((Q-1):-1:1)]/2);
    L = length(X);
    
    for i = 1:S(1)
        if mod(i,round(S(1)/4))==0; disp(num2str(i)); end
        dum = input(i,:,:);
        parfor(j = 1:S(3), PoolSize)
            input(i,:,j) = HaarTransformIntMultiscale_v1p1(dum(1,:,j),X,L,pad);
        end
    end
elseif dim==2
    
    N = numel(input(:,1,1));
    pad = mod(N,2);
    Q = nextpow2(N);
    X = round([N 2.^((Q-1):-1:1)]/2);
    L = length(X);
    
    for i=1:S(2)
        if mod(i,round(S(2)/4))==0; disp(num2str(i)); end
        dum = input(:,i,:);
        parfor(j = 1:S(3), PoolSize)
            input(:,i,j) = HaarTransformIntMultiscale_v1p1(dum(:,1,j),X,L,pad);
        end
    end
elseif dim==3
    disp('Haar on 3rd dimension not tested!  And it''s usually slow!')
    
    N = numel(input(1,1,:));
    pad = mod(N,2);
    Q = nextpow2(N);
    X = round([N 2.^((Q-1):-1:1)]/2);
    L = length(X);
    
    for i=1:S(1)
        dum = input(i,:,:);        
        parfor(j = 1:S(2), PoolSize)
            input(i,j,:) = HaarTransformIntMultiscale_v1p1(dum(1,j,:),X,L,pad);
        end
    end
end



function input = HaarTransformIntMultiscale_v1p1(input,X,L,pad)
%% Comments
% input = HaarTransformIntMultiscale_v1p1(input,X,L,pad)
%
% This function takes in an 'input' vector and returns it after performing
% a fully iterated Haar transform on successively smaller scales (i.e.,
% subsampled 'blur' sections) until it reaches the smallest possible size.
% It accepts arguments:
%   input:  the vector to be transformed
%   X:      scales at which to perform the transform
%   L:      the number of such scales
%   pad:    binary switch; 1 = pad w/ zeros to next power of 2, 2 = don't.

if nargin==1
    clc
    N = numel(input);
    pad = mod(N,2);
    Q = nextpow2(N);
    X = round([N 2.^((Q-1):-1:1)]/2);
    L = length(X);
    input;
end

%% Do Haar transform over successively smaller pieces
% If input has an odd number of elements, pad with an extra zero
i0 = 1;
if pad
    input(end+1) = 0;
    
    input = reshape(input,2,X(1));
    input(2,:) = input(2,:) - input(1,:);
    input(1,1:end-1) = input(1,1:end-1) + input(2,1:end-1)/2;
    i0 = 2;
    input = input';
    input = input(:);
end

for i=i0:L
    input(1:2*X(i)) = HaarTransformInt_v1p1(input(1:2*X(i)),2*X(i));
end

if pad
    input(end) = [];
end

function input = HaarTransformInt_v1p1(input,K)
%%
% input = HaarTransformIntMultiscale(input,K)
%
% This function takes in an 'input' vector and returns it after performing
% a single Haar transform on the largest possible scale (i.e., the whole
% vector).
% It accepts arguments:
%   input:  the vector to be transformed
%   K:      the length of the vector
%
% The reasons for passing in the length is
% because this function is often called as a subroutine for large datasets
% where it may be called millions of times.  Thus, it pays to keep
% operations within this function to a minimum.

% Reshape into 2 rows and length/2 columns to facilitate Haar transform
input=reshape(input,2,K/2);
% input = [input(1:2:end)'; input(2:2:end)'];

% dum = zeros(2,K/2);
% dum(:) = input(:);
% input = dum;


% Use lifting technique from Calderbank 1998 to maintain integer datatype
% integrity during this operation
input(2,:) = input(2,:) - input(1,:);
input(1,:) = input(1,:) + input(2,:)/2;

% Reshape vector to return in correct column format
input = input';
input = input(:);


