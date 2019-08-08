    function PlotCirc(img,a,b,R,h,color)
    
    %% Comments
    % This function takes in an image with a circle at center (a,b) with
    % radius R and plots this image with a fitted circle in figure h.
    
    %%
%     if h==gcf
%         1
%         h = gca
%     end
    t = 0:0.05:2*pi+.05;
       x = a + R*cos(t);
       y = b + R*sin(t);
       pcolor(h,img);
       daspect(h,[1 1 1]);
       shading(h,'flat');
%        hold(h);
       line(x,y,'Parent',h,'Color',color,'LineWidth',3)
       line(a,b,'Parent',h,'Color','y','Marker','.','MarkerSize',16)
%        hold(h);
    end