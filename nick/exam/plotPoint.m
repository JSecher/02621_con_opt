function [plt] = plotPoint(x, y, type, color, sz)
%plotPoint Plots a point on a countur plot
%   Type gives the marker

if nargin < 3
    type = "gen"; % generic
end

if type == "max" % Maximum
    a = '^'; b='MarkerFaceColor'; c='r'; d='MarkerEdgeColor'; e='k'; f='markersize';g=15;
elseif type == "min" % Minimum
    a = 'v'; b='MarkerFaceColor'; c='g'; d='MarkerEdgeColor'; e='k'; f='markersize';g=15;
elseif type == "sad" % Saddle
    a = 'd'; b='MarkerFaceColor'; c='m'; d='MarkerEdgeColor'; e='k'; f='markersize';g=15;
elseif type == "gen" % Just some point
    a = 'o'; b='MarkerFaceColor'; c='b'; d='MarkerEdgeColor'; e='k'; f='markersize';g=15;
elseif type == "int" % A inital starting point
    a = 's'; b='MarkerFaceColor'; c='c'; d='MarkerEdgeColor'; e='k'; f='markersize';g=15;
elseif type == "sol" % Found solution point
    a = 'h'; b='MarkerFaceColor'; c='y'; d='MarkerEdgeColor'; e='k'; f='markersize';g=15;
else

    error("Unknown type of point: %s", type)
end

if nargin >= 4
    c = color;
end


if nargin >= 5
    g = sz;
end

plt = plot(x,y,a,b,c,d,e,f,g);

end