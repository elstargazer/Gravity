function figure1=PlotHPowers(Y1)
%CREATEFIGURE(Y1)
%  Y1:  vector of y data

%  Auto-generated by MATLAB on 05-Jul-2012 20:27:19

% Create figure
figure1 = figure('XVisual','','Color',[0 0 0]);

% Create axes
axes1 = axes('Parent',figure1,'ZColor',[1 1 1],'YScale','log',...
    'YMinorTick','on',...
    'YColor',[1 1 1],...
    'XColor',[1 1 1],...
    'FontSize',30,...
    'FontName','Times New Roman',...
    'Color',[0 0 0]);
box(axes1,'on');
hold(axes1,'all');

% Create semilogy
semilogy(Y1,'MarkerFaceColor',[0 0 1],'MarkerSize',4,'Marker','o',...
    'LineWidth',2);

grid on;


ylabel('Total power');
xlabel('Degree');
