%% https://dgleich.github.io/hq-matlab-figs/
% The new defaults will not take effect if there are any open figures. To
% use them, we close all figures, and then repeat the first example.
close all;

width = 4.74;%3.19;%2.96;%5.93;%3.3;     % Width in inches
height = 4.74;%2.20;%3.3;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 9;%11;      % Fontsize
lw = 0.5; %1.5;      % LineWidth
msz = 3;%8;       % MarkerSize

% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
set(0,'defaultErrorBarLineWidth',lw);   % set the default line width to lw
set(0,'defaultErrorBarMarkerSize',msz); % set the default line marker size to msz
%set(0,'defaultAxesFontSize',fsz);

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
set(0, 'defaultFigurePaperSize', [width+0.2 height+0.2]);
set(0, 'defaultFigurePaperPosition', [0.1 0.1 width height]);

set(0, 'defaultTextInterpreter', 'latex');
set(0, 'defaultFigureRenderer', 'painters');

%% Save
print('pure_tone','-dpdf','-r300');
