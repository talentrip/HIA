function dispmat = SpokeSimulator(ImagePixels,ChannelPixelRadius,ChannelPixelWidth,angfreq,SpokeWidthDeg,NumFrames,EdgeBrightnessRatio)

%% Comments
% function gerbilsimv3p8(ChannelPixelRadius,ChannelPixelWidth,angfreq,SpokeWidthDeg,framrat,NumFrames,EBR)
% simulates a gerbil in a Hall Ion Thruster in movie form
%
% Inputs:
%
% ChannelPixelRadius = mean radius in pixels, total size is 160 x 160
% ChannelPixelWidth = mean channel diameter in pixels
% angfreq = angular frequency in degrees theta per frame
% SpokeWidthDeg = gerbil size in degrees theta
% framrat = frame rate in fps for eventual movie processing
% NumFrames = total desired movie time in seconds
% EBR = fractional edge brightness ratio for use with gaussian gerbil profile

%sample input gerbilsimv3p7(45,15,4,14,40,20,.15)


%% Constants and switches

%tic
MaxNumSpokes = 3;

%theta0 will be the vector containing the gerbils theta position
%NOTE the length of theta0 will tell us how many gerbils are present at a
%specific time (NumFrames), if the position is NaN that # gerbil is not
%present
%no gerbils present at start of simulation
CurrentSpokeAngle  = NaN(1,MaxNumSpokes);
thetalo = NaN(1,MaxNumSpokes);
thetahi = NaN(1,MaxNumSpokes);

%gerblongev: number of frames first gerbil will loop through
%If first gerbil is to loop through 50 times this corresponds to a theta
%longevity of 50*angfreq, if the gerblongev is NaN then that position
%gerbil is not alive

CurrentNumSpokes = 0;                                                             % start num of gerbils; later, current num of gerbils

SpokeLifeStdDev = 60;    %standard deviation used for longevity calculation (number of frames)
SpokeLifeMean = 360;   %mean theta distance used in longevity calculation

gaussprofile = 1;                                                          % set to 0 for step-fn borders
debug = 0;
debug2 = 0;

%% Initialize other variables

dispmat = zeros(ImagePixels,ImagePixels,NumFrames);                                          %preallocates matrices needed for final gerbil display
gaussian = zeros(ImagePixels,ImagePixels,3);
thetapicker = zeros(ImagePixels,ImagePixels,3);

x0 = ImagePixels/2;
y0 = ImagePixels/2; 
[r theta] = PixelPolarCoords(dispmat(:,:,1),[x0 y0]);

ChannelIR = ChannelPixelRadius - ChannelPixelWidth/2;
ChannelOR = ChannelPixelRadius + ChannelPixelWidth/2;
%makes a logical matrix of a ring (ring part is 1, non ring 0)
ring = r>ChannelIR & r<ChannelOR;

%% determine decay constant k for Gaussian

if(gaussprofile)
   % EdgeBrightnessRatio = 0.1;
   k = -log(EdgeBrightnessRatio)*(2/SpokeWidthDeg)^2;
else
   k = 0;
end

%% create predetermined gerbil vectors

SpokeLifeDist = randn(1,NumFrames)*SpokeLifeStdDev + SpokeLifeMean;                                       % vector w/ preditermined theta longevity values based on a normal distribution
SpokeLifeDist = round(SpokeLifeDist);                                                % round to integer value to use as a number of loop iterations to survive
SpokeLifetimes = SpokeLifeDist./angfreq;                                                % vector of predetermined NumFramess per gerbil
% randspot = 360*rand(1,NumFrames) - 180;                                      % randspot is a vector containing 5000 random numbers between -180 & 180
CurrentSpokeLifetime = SpokeLifetimes(1:3);
%vector of predetermined starting
%locations for gerbil number one

   gerb1prob = 0.15;
   gerb2prob = 0.35;
   gerb3prob = 0.55;

%% loop through movie frames, creating gerbil images

counter = [0 0 0]; %initialize counter
cumdeaths = 0;
for i = 1:NumFrames

   %% determine if a gerbil should die or not
   counter = counter+1;
   death = counter > CurrentSpokeLifetime;
   counter(death) = 0;

   numdeaths = sum(death);
   cumdeaths = cumdeaths + numdeaths;
   CurrentNumSpokes = CurrentNumSpokes - numdeaths;

   CurrentSpokeAngle(death) = NaN;
   dead = find(isnan(CurrentSpokeAngle));
   live = find(~isnan(CurrentSpokeAngle));

   CurrentSpokeLifetime(death==1) = SpokeLifetimes(cumdeaths+MaxNumSpokes+1:cumdeaths+MaxNumSpokes+numdeaths);    

   while (max(CurrentSpokeAngle) > 180)                                 % This while statement occurs often in the code --> It makes sure
       CurrentSpokeAngle(CurrentSpokeAngle>180) = CurrentSpokeAngle(CurrentSpokeAngle>180)-360;                     % the theta value assigned is inbetween -180 & 180
   end


   %% probability scenario

   %The generation of a gerbil will be determined by a probability scenarion.
   %The probability that a gerbil will generate depends on how many gerbils
   %are currently present.  The location of the new gerbil will depend on
   %where the other gerbils are present.  A new gerbil's position will be
   %located the farthest distance from all other gerbils.

   %Note to user: If you wish to increase the probability of a gerbil being created, then
   %decrease the variable 'gerb#prob' and visa versa

   switch CurrentNumSpokes
       
       case 0   % Called if no gerbils are currently present
       if (rand >= gerb1prob)                                         % if true create new gerbil
           CurrentSpokeAngle(1) = 360*rand(1)-180;                                   % new gerbil is created at a random theta location
           thetalo(1) = 0;     thetahi(1) = 0;
           counter(1) = 0;
           CurrentNumSpokes = CurrentNumSpokes + 1;
           CurrentSpokeLifetime(1) = SpokeLifetimes(i);                                 % selects a longivitey(number of frames gerbil will live) for the corresponding gerbil
       end

       case 1   % Called if one gerbil is currently present
       if (rand >= gerb2prob)                                        % if true create 2nd gerbil diametrically opposite to first
           if (debug2)
%                 [live, dead] = indexgerbalive(CurrentSpokeAngle);                     % indexgerbalive is a called sub function that determines the index of the live and dead gerbil
               live
               dead
           end
           CurrentSpokeAngle(dead(1)) = CurrentSpokeAngle(live) + 180;                         % the location of this new gerbil will be on excact opposite side of current gerbil
           thetalo(dead(1)) = 0;      thetahi(dead(1)) = 0;
           CurrentNumSpokes = CurrentNumSpokes + 1;
           counter(dead(1)) = 0;
           CurrentSpokeLifetime(dead(1)) = SpokeLifetimes(i);                              % selects a longivitey(number of frames gerbil will live) for the corresponding gerbil
       end

       case 2                                              % Called if two gerbils are currently present
       if (rand >= gerb3prob)                                         % if true create 3rd gerbil as far as possible from other two
           %disp('INVOKING THIRD GERBIL PLACMENT ALGORITHM')
           livethetas = CurrentSpokeAngle(live);                               % determines the theta values of the live gerbils
           spot = [0 0];                                              % next several lines ensure new gerbil forms as far as possible from old ones
           spot(1) = mean(livethetas);
           if spot(1)>0
               spot(2) = spot(1)-180;
           else
               spot(2) = spot(1)+180;
           end      
           gaps = livethetas(1)-spot;
           if(max(abs(gaps))>180)
               gaps(gaps>180) = 360 - gaps(gaps>180);
               gaps(gaps<-180) = 360 + gaps(gaps<-180);
           end
           [gap,which] = max(abs(gaps));
           CurrentSpokeAngle(dead) = spot(which);
           CurrentNumSpokes = CurrentNumSpokes+1;            
       end
       
   end       %end of CurrentNumSpokes if statement

   %% contiuation of main for loop (looping for through num frames, creating gerbil image)

   CurrentSpokeAngle = CurrentSpokeAngle + angfreq;             %sets a new theta position to all gerbils

   for v = 1:3
       thetalo(v) = CurrentSpokeAngle(v) - SpokeWidthDeg/2;             %calculates lower bound of gerbil
       while (thetalo(v) > 180)
           thetalo(v) = thetalo(v) - 360;           %makes sure theta is not out of the bounds for our circle
       end
       thetahi(v) = thetalo(v) + SpokeWidthDeg;              %calculates upper bound of gerbil
       while (thetahi(v) > 180)
           thetahi(v) = thetahi(v) - 360;
       end
   end


   %% note from Mike: try to generalize logic here, perhaps by looking for max values of
   % the sum of the two halves of thetapicker
   for z = 1:3     %loop through number of gerbils times

       ifnan = isnan(CurrentSpokeAngle);       %determines what gerbils are designated as a NaN

       if (~ifnan(z))   %called if a gerbil is present at that theta index

           %compentation for the jump.  Occurs when thetalo is a positive number in
           %the second quadrant and thetahi is a negative number in the third
           %quadrant

           if (((thetalo(z) >= 90)&&(thetalo(z) <= 180))&&((thetahi(z) >= -180)&&(thetahi(z) <= 90)))%called if thetalo is in 2nd quadrant &
               %thetahi is in 3rd quadrant
               thetapickerlo(:,:,z) = (theta > thetalo(z));
               thetapickerhi(:,:,z) = (theta < thetahi(z));
               thetapicker(:,:,z) = thetapickerlo(:,:,z) + thetapickerhi(:,:,z);%add thetapickers to compensate for jump

               gaussianlo(:,:,z) = exp ( -k* (theta - (thetahi(z) + thetalo(z) + 360)/2).^2);
               gaussianhi(:,:,z) = exp ( -k* (theta - (thetahi(z) + thetalo(z) - 360)/2).^2);
               gaussian(:,:,z) = gaussianlo(:,:,z) + gaussianhi(:,:,z);         %add gaussian hi and lo to compensate for jump

           else                                                                                   %called if the scenario does NOT call
               %for a compensation
               thetapicker(:,:,z) = (theta>thetalo(z)).*(theta<thetahi(z));%sets a logical matrix equal to just the gerbil
               gaussian(:,:,z) = exp ( -k* (theta - (thetahi(z) + thetalo(z))/2).^2);%applies gaussian function to gerbil

           end
       else
           thetapicker(:,:,z) = zeros(ImagePixels,ImagePixels);     %the thetapicker and gaussian corresponding to a dead gerbil                                        %will just be a matrix of all zeros
           gaussian(:,:,z) = zeros(ImagePixels,ImagePixels);        %will be a matrix of all zeros
       end

   end  %ends nested for loop 1:3


   if(debug)
       [thetalo thetahi]

       delay = .25;
       pcolor(ring)
       colormap gray;
       shading interp
       daspect([1 1 1])

       pause(delay)
       pcolor(thetapicker)
       colormap gray;
       shading interp
       daspect([1 1 1])
       pause
       pcolor(gaussian)
       colormap gray;
       shading interp
       daspect([1 1 1])
       pause
   end


   %the following will layer all of the image matrices for each gerbil
   %ontop of each other.  gaussians and thetapickers are the final matrices
   %that are used to make the output matrix for this particular NumFrames.
   %gaussians and thetapickers depend on the number of gerbils present
   %aka CurrentNumSpokes

   gaussians = sum(gaussian,3);
   thetapickers = sum(thetapicker,3);

   dispmat(:,:,i) = ring.*gaussians.*thetapickers;  %connect the puzzle pieces


%    % plotting options
%    pcolor(dispmat(:,:,i))%displays the gerbil at location thetapicker(:,:,i)
%    colormap gray;
%    shading interp        %make gerbil look prettier
%    daspect([1 1 1])
%    pause(.1/framrat)     %helps determine the speed of gerbil

end  %end to big NumFrames loop
%toc

%% normalize dispmat to 8-bit color for movie frame conversion
% find peak brightness over all bitmaps
high = max(max(max(abs(dispmat))));
% rescale to 8-bit (256 color) resolution
dispmat = dispmat /(high+1) * 2^5;
% make sure all elements are >0 (Zero elements give errors in im2frame)
dispmat = dispmat + 2^5 + 1;

% %% make the movie
% for i = NumFrames:-1:1
%    M(i) = im2frame(dispmat(:,:,i), gray);    %layers on all frames to make a movie
% end
% figure(2)
% movie(M)        %shows newly created movie
% 
% movie2avi(M,'gerbilmovie.avi');     %creates a .avi movie file of simulation, will be created in current directory with name gerbilmovie.avi