function [c, con] = contourHimmel(add_con, view, xminmax, yminmax, res, f, varargin)
% contourplot Generates a contour plot of the given function
% 
% Output:
%   figure    : The figure object
%%

if nargin < 5
    f = @(x1,x2) (x1.^2+x2-11).^2 + (x1+x2.^2-7).^2;
end

if nargin < 1
    add_con = false;
end

if nargin < 2
    view = 6;
end

if nargin < 3
    xminmax = 5;
    yminmax = 5;
    res = 50;
end

x = -view:1/res:view;
y = -view:1/res:view;
[X,Y] = meshgrid(x,y);
v = [0:2:10 10:10:100 100:20:300];
F = f(X,Y);
[m,c]=contour(X,Y,F,v,"linewidth",2);
colormap("turbo")

% Add constraints 
con = [];
if add_con
    % Con 1
    p1 = (x+2).^2;
    idx = abs(p1) <= view;
    p1 = p1(idx);
    xp = x(idx);
    xp = [xp,xp(end),xp(1),xp(1)];
    p1 = [p1,view, view, p1(1)];
    con1 = patch("XData", xp,"YData", p1,'FaceColor', 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    
    % Con 2
    p2 = 4/10*x;
    idx2 = abs(x) <= view;
    p2 = p2(idx2);
    xp2 = x(idx2);
    xp2 = [xp2,xp2(end),xp2(1),xp2(1)];
    p2 = [p2,-view,-view,xp2(1)];
    con2 = patch("XData", xp2,"YData", p2,'FaceColor', 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    
    % X con
    con3 = patch("XData", [-view, -xminmax, -xminmax,-view, -view],...
          "YData", [-view, -view, view, view, -view],...
          'FaceColor', 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    con4 = patch("XData", [xminmax, view, view, xminmax, xminmax],...
          "YData", [-view, -view, view, view, -view],...
          'FaceColor', 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k');

    % Y con
    con5 = patch("XData", [-view, -view, view, view, -view],...
          "YData", [view, yminmax, yminmax, view, view],...
          'FaceColor', 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
    con6 = patch("XData", [-view, -view, view, view, -view],...
            "YData", [-yminmax, -view, -view, -yminmax, -yminmax],...
            'FaceColor', 'red', 'FaceAlpha', 0.5, 'EdgeColor', 'k');

    con = [con1, con2, con3, con4, con5, con6];

end
        
% Turns off the legend for the contour
set(get(get(c,'Annotation'),'LegendInformation' ),'IconDisplayStyle', 'off' );

% Appropiate labels for the plot
xlabel('x_1','Fontsize',14)
ylabel('x_2','Fontsize',14)
colorbar

axis([-view view -view view])

end