function [RMSvalue] = rms(vector,dim)
%Compute the Root-Mean-Square of a vector

if nargin==2 && dim==2
    RMSvalue    = sqrt(mean(vector.^2,2));
else  
    RMSvalue    = sqrt(mean(vector.^2));    
end
