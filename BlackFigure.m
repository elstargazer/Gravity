function [figure1,axes1]=BlackFigure(ttl,ylbl,xlbl)
%CREATEFIGURE(X1,Y1)
%  X1:  vector of x data
%  Y1:  vector of y data

%  Auto-generated by MATLAB on 17-Oct-2012 17:12:26

% Create figure
figure1 = figure('XVisual','','Color',[0 0 0]);

% Create axes
axes1 = axes('Parent',figure1,'ZColor',[1 1 1],...
    'YMinorTick','on',...
    'YMinorGrid','on',...
    'YColor',[1 1 1],...
    'XColor',[1 1 1],...
    'FontSize',20,...
    'FontName','times',...
    'Color',[0 0 0]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');

% Create ylabel
ylabel(ylbl,'FontSize',25,'FontName','times',...
    'Color',[1 1 1]);

% Create xlabel
xlabel(xlbl,'FontSize',25,'FontName','times','Color',[1 1 1]);

% Create title
title(ttl,'FontSize',25,'FontName','times',...
    'Color',[1 1 1]);

