function [c, con] = contourRosen(add_con, xview,yview, xminmax, yminmax, res, f, varargin)
% contourRosen Generates a contour plot of the given function
% 
% Output:
%   figure    : The figure object
%%

if nargin < 5
    f = @(x1,x2) 100 * (x2-x1.^2).^2 + (1-x1).^2;
end

if nargin < 1
    add_con = 0;
end

if nargin < 2
    xviewmin = -1.6;
    xviewmax = 1.6;
else
    xviewmin = min(xview);
    xviewmax = max(xview);
end

if nargin < 3
    yviewmin = -1;
    yviewmax = 2;
else
    yviewmin = min(yview);
    yviewmax = max(yview);
end

absviewmax = max(abs([xviewmin, xviewmax, yviewmin, yviewmax]));

if nargin < 4
    xminmax = 2;
    yminmax = 2;
    res = 50;
end

x = xviewmin:1/res:xviewmax;
y = yviewmin:1/res:yviewmax;
[X,Y] = meshgrid(x,y);
v = [0:0.1:2 2:1:10 10:10:300];
F = f(X,Y);
[m,c]=contour(X,Y,F,v,"linewidth",2);
colormap("turbo")

% Add constraints 
con = [];
if add_con == 1
    % Circular constraint

    % Get data for the two circles
    [xp5, yp5] = circle(0,0,1);
    [xp10,yp10] = circle(0,0,10);
    % Draw the two circles
    hold on
    plot(xp5,yp5)
    plot(xp10,yp10)
    % Fill in the area by defining the polygon that is defined by tracing the inner circle and then back along the outer circle
    h = fill([xp5 flip(xp10)],[yp5 flip(yp10)],'g');
    set(h,'facealpha',.5,'FaceColor', 'red', 'EdgeColor', "k")
    % Cover up the connecting line in the same color
    h1 = line([xp5(1) xp10(end)],[yp5(1) yp10(end)]);
    set(h1,'Color','r')
    % Set the axes to have equal length, so that they *look* like circles
    axis equal
    con = [h];

elseif add_con == 2
    % box constraint
    xb = [-1.5; 1.5];      % Lower bound for x
    yb = [0.5; 1.5];        % Upper bound for x
    

    % X con
    con3 = patch("XData", [-absviewmax, xb(1), xb(1),-absviewmax, -absviewmax],...
          "YData", [-absviewmax, -absviewmax, absviewmax, absviewmax, -absviewmax],...
          'FaceColor', 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    con4 = patch("XData", [xb(2), absviewmax, absviewmax, xb(2), xb(2)],...
          "YData", [-absviewmax, -absviewmax, absviewmax, absviewmax, -absviewmax],...
          'FaceColor', 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k');

    % Y con
    con5 = patch("XData", [-absviewmax, -absviewmax, absviewmax, absviewmax, -absviewmax],...
          "YData", [absviewmax, yb(2), yb(2), absviewmax, absviewmax],...
          'FaceColor', 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    con6 = patch("XData", [-absviewmax, -absviewmax, absviewmax, absviewmax, -absviewmax],...
            "YData", [yb(1), -absviewmax, -absviewmax, yb(1), yb(1)],...
            'FaceColor', 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    axis equal
    con = [con3, con4, con5, con6];

elseif add_con == 3
    % box constraint
    xb = [-1.5; 1.5];      % Lower bound for x
    yb = [0.5; 1.5];        % Upper bound for x
    

    % X con
    con3 = patch("XData", [-absviewmax, xb(1), xb(1),-absviewmax, -absviewmax],...
          "YData", [-absviewmax, -absviewmax, absviewmax, absviewmax, -absviewmax],...
          'FaceColor', 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    con4 = patch("XData", [xb(2), absviewmax, absviewmax, xb(2), xb(2)],...
          "YData", [-absviewmax, -absviewmax, absviewmax, absviewmax, -absviewmax],...
          'FaceColor', 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k');

    % Y con
    con5 = patch("XData", [-absviewmax, -absviewmax, absviewmax, absviewmax, -absviewmax],...
          "YData", [absviewmax, yb(2), yb(2), absviewmax, absviewmax],...
          'FaceColor', 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    con6 = patch("XData", [-absviewmax, -absviewmax, absviewmax, absviewmax, -absviewmax],...
            "YData", [yb(1), -absviewmax, -absviewmax, yb(1), yb(1)],...
            'FaceColor', 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    axis equal
    con = [con3, con4, con5, con6];

        % Con 1
    p1 = -(2*x).^2+1;
    idx = abs(p1) >= yviewmin;
    p1 = p1(idx);
    xp = x(idx);
    xp = [xp,xp(end),xp(1),xp(1)];
    p1 = [p1,yviewmin, yviewmin, p1(1)];
    con1 = patch("XData", xp,"YData", p1,'FaceColor', 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    
    con = [con3, con4, con5, con6,con1];
end
        
% Turns off the legend for the contour
set(get(get(c,'Annotation'),'LegendInformation' ),'IconDisplayStyle', 'off' );

% Appropiate labels for the plot
xlabel('x_1','Fontsize',14)
ylabel('x_2','Fontsize',14)
colorbar

axis([xviewmin xviewmax yviewmin yviewmax])

end


function [xunit, yunit] = circle(x,y,r)
    
    th = 0:pi/50:2*pi;
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
end
    