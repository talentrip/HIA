function Data = VideoDecompress_v1p1(Data)

% Order of operations:  Haar compression of rows and columns will normally
% commute, i.e., you can transform across rows, then columns, rows,
% columns, etc. and inverse transform back again without keeping track of
% order.
%
% However, this only works in floating precision (Single or double)
% arithmetic.  In integer arithmetic, order matters and you can introduce
% roundoff errors all over the place if you're not careful.  For this
% reason, the steps below perform column inverse transforms first, followed
% by rows, in order to match with the order used in the video compression
% (which was rows first, then columns).

% Optional steps for images
imageopt = 0;
if imageopt
    MeanData = mean(Data,3,ParallelOK);
    % normalization?
    
    Data = bsxfun(@minus,Data,MeanData);
end
tic


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

% Haar transform columns
disp('Undoing column Haar...')
Data = HaarInverseMultiDim_v1p1(Data,2,ParallelOK);
toc

% Rows
disp('Undoing row Haar...')
Data = HaarInverseMultiDim_v1p1(Data,1,ParallelOK);
toc

% Stacks: Don't do stacks
% if size(Data,3)~=1
%     disp([sprintf('\r') 'Undoing 3D Haar...' sprintf('\r')])
%     Data = HaarInverseMultiDim(Data,3);
% end
% toc

function input = HaarInverseMultiDim_v1p1(input,dim,ParallelOK)

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
    
    for i=1:S(1)
        if mod(i,round(S(1)/10))==0; disp(num2str(i)); end
        dum = input(i,:,:);
        parfor(j = 1:S(3), PoolSize)
            input(i,:,j) = HaarTransformInverseIntMultiscale_v1p1(dum(1,:,j),X,L,pad);
        end
    end
elseif dim==2
    
    N = numel(input(:,1,1));
    pad = mod(N,2);
    Q = nextpow2(N);
    X = round([N 2.^((Q-1):-1:1)]/2);
    L = length(X);
    
    for i=1:S(2)
        if mod(i,round(S(2)/10))==0; disp(num2str(i)); end
        dum = input(:,i,:);
        parfor(j = 1:S(3), PoolSize)
            input(:,i,j) = HaarTransformInverseIntMultiscale_v1p1(dum(:,1,j),X,L,pad);
        end
    end
elseif dim==3
    
    N = numel(input(1,1,:));
    pad = mod(N,2);
    Q = nextpow2(N);
    X = round([N 2.^((Q-1):-1:1)]/2);
    L = length(X);
    
    for i=1:S(1)
        for j = 1:S(2)
            input(i,j,:) = HaarTransformInverseIntMultiscale_v1p1(input(i,j,:),X,L,pad);
        end
    end
end

Sparse8BitFraction = sum(abs(input(:))>2^7)/numel(input);
disp([num2str(100*Sparse8BitFraction) '% of image pixels require > 8 bits' sprintf('\r')])

function input = HaarTransformInverseIntMultiscale_v1p1(input,X,L,pad)
%% Comments
% input = HaarTransformIntMultiscale(input,X,L,pad)
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

%% Do the successively smaller Haar transforms
% If input has an odd number of elements, pad with an extra zero
i0 = 1;
if pad
    input(end+1) = 0;
end

for i=L:-1:i0
    input(1:2*X(i)) = HaarTransformInverseInt_v1p1(input(1:2*X(i)),2*X(i));
end

if pad
    input(end) = [];
end

function input = HaarTransformInverseInt_v1p1(input,K)
% if length(size(input))>2 || min(size(input))>1
%    disp('Error, cannot take multidimensional inputs yet')
%    return
% end

input = reshape(input,K/2,2);

% "Lifting" technique
input(:,1) = input(:,1) - input(:,2)/2;
input(:,2) = input(:,1) + input(:,2);

input = input';
input = input(:);






