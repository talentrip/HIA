function isaSubroutine = CheckIfSubroutine(thisFile)

% Pull up the hierarchy of all calling functions using 'dbstack'
[a,b] = dbstack;
% The file string at the end of 'a.name' is the topmost parent calling
% function
callingFile = a(end).name;

% If that topmost parent calling function is the same one that queried this
% script, then it is not a subroutine.  Otherwise, it is a subroutine.
if strcmp(callingFile,thisFile)
    isaSubroutine = 0;
else
    isaSubroutine = 1;
end