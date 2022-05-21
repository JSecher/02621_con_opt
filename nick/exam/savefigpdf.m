function [path] = savefigpdf(fig, filename, ex)
%savefigpdf Save figure to a pdf of right size
%   
% Inputs:
%   fig         : Figure object to save
%   filename    : Filename to be saved to, without extension
%   ex          : Excersice number, integer from 1 to 5, giving the folder to
%                 save to
%
% Outputs:
%   path        : a string giving the path to where the file was saved

foldername = "../graphics";
if ~exist(foldername, 'dir')
       mkdir(foldername);
end
if nargin > 2
    foldername = sprintf("%s/ex%d/", foldername, ex);
    if ~exist(foldername, 'dir')
           mkdir(foldername);
    end
end

path = foldername + filename;
pos = get(fig,'Position');
set(fig,'PaperPositionMode','Auto','PaperUnits','points','PaperSize',[pos(3), pos(4)]);
set(gca,'FontSize',16) 
print(fig, path, '-dpdf','-fillpage');
fprintf('... Plot has been saved to ''%s''\n', path);
end

