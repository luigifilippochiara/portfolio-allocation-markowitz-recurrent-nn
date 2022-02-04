% source https://it.mathworks.com/matlabcentral/fileexchange/34087-cbrewer-colorbrewer-schemes-for-matlab
% useful coloring techniques: http://math.loyola.edu/~loberbro/matlab/html/colorsInMatlab.html

% set default color palette for plots
% 20=nr of assets (it's mandatory to set the number of colors)
% C = cbrewer('qual','Paired',20);
C = cbrewer('qual','Set1',20);
set(0,'DefaultAxesColorOrder',C)
% set default line width
set(0,'DefaultLineLineWidth',1.2)
set(0,'DefaultFigureColormap',cbrewer('seq','YlOrRd',64));
set(0,'DefaultAxesTitleFontSizeMultiplier',2)
set(0,'DefaultAxesFontSize',13)