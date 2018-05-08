% create some data to plot
t = 1:100;
y1 = abs(1e4*sin(.01*t));
y2 = abs(cos(.02*t));
%% create a figure and use PLOTYY, saving the output handles
fg = figure(5);
% note that by default the tick marks overlap exactly (and it seems that
% only one set of tick marks is present per line axis
[AX,H1,H2] = plotyy(t, y1, t, y2);
% press ENTER to continue
pause;
%% change the number of tick marks on each Y axis to see the issue
Nticks = 8;
y1 = linspace(0, 8e4, Nticks);
y2 = linspace(0, 14, Nticks+4); %put 4 more ticks on the right Y axis
% set the new limits and tick marks. Notice how there appear to be two sets
% of tick marks on the right axis.
set(AX(1), 'ylim', [y1(1), y1(end)], 'ytick', y1);
set(AX(2), 'ylim', [y2(1), y2(end)], 'ytick', y2); 
% press ENTER to continue
pause;
%% link the two X axes
linkaxes(AX,'x');
%% turn off the axis Box property for box axis
% the box property results in the display of tick marks on both axis
set(AX(1),'Box','off')
set(AX(2),'Box','off')
% notice how the bounding box on top of the image is missing
%% replace the top edge of the box with the second X axes
% and turn off its labels
set(AX(2), 'XTickLabel','','XAxisLocation','Top') 
%% print the figure to .png format
% to confirm that the output is as desired
print(fg,'-dpng','myFigure');