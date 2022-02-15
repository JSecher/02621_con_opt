%% Problem 2 - Linear Optimization


function ccontourplot(f, ...
                      xmin, xmax, ymin, ymax, ...
                      xres, yres, ...
                      varargin)

    x = xmin:xres:xmax;
    y = ymin:yres:ymax;
    [X,Y] = meshgrid(x,y);
    
    v = [0:2:10 10:10:100 100:20:200]
    [c,h]=contour(X,Y,f,v,"flinewidth",2);
    colorbar
    
    yc1 = (x+2).^2;
    yc2 = (4*x)/10;
    hold on
        fill(x,yc1,[0.7 0.7 0.7],’facealpha’,0.2)
        fill([x x(end) x(1)],[yc2 -5 -5],[0.7 0.7 0.7],’facealpha’,0.2)
    hold off

end

