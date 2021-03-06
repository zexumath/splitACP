function wysiwyg
%WYSIWYG -- this function is called with no args and merely
%%       changes the size of the figure on the screen to equal
%%       the size of the figure that would be printed, 
%%       according to the papersize attribute.  Use this function
%%       to give a more accurate picture of what will be 
%%       printed.
%
%%   Dan(K) Braithwaite, Dept. of Hydrology U.of.A  11/93
%%   Modified by Toby Driscoll, University of Delaware, 2003.  
% 

% set(gcf, 'PaperPositionMode', 'manual');
% set( gcf, 'PaperUnits', 'centimeters' );
% set( gcf, 'PaperType', 'usletter' );
% width = 12.0;
% height = 12.0;
% set( gcf, 'PaperPosition', [0, 0, width, height] );
% set(gca, 'FontSize', 14);
% set(gca, 'FontName', 'Times New Roman', 'FontWeight','Bold');
% box on

unis = get(gcf,'units');
ppos = get(gcf,'paperposition');
set(gcf,'units',get(gcf,'paperunits'));
pos = get(gcf,'position');
pos(2) = max(0,pos(2)-1.05*(ppos(4)-pos(4)));
pos(3:4) = ppos(3:4);
set(gcf,'position',pos);
set(gcf,'units',unis);

