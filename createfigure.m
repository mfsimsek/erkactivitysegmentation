%open the lut for SFC colors
load('matlab_SFCLUT.mat');

%rescale SFC values to 0 - 0.44 range min to max in special LUT
SFC0=SFC/0.44;
SFC0(SFC0>1)=1;
SFC1=SFC0;
SFC1(isnan(SFC1))=0;
H = fspecial('disk',3);
SFC1 = imfilter(SFC1,H,'replicate'); 
SFC1(isnan(SFC0))=NaN;
figure
mesh(SFC1(:,T-380:T))

%use below to convert any mesh plot into 2-D figures
fig = gcf;
axObjs = gca;
dataObjs = axObjs.Children;
x = dataObjs(1).XData;
y = dataObjs(1).YData;
zdata1 = dataObjs(1).ZData;

% Create figure
figure('InvertHardcopy','off','Color',[1 1 1],...
    'OuterPosition',[276 287 576 513]);
colormap(LUTSFC); %import the lut first!

% Create axes
axes1=gca;

% Create surf
surf1 = surf(zdata1,zdata1,'FaceLighting','none','EdgeLighting','flat',...
    'FaceColor',[1 1 1],...
    'EdgeColor','interp');


% Create ylabel
ylabel('PSM (cells)');

% Create xlabel
xlabel('Time (min)');

% Create title
%%use for different inhibitory drug pulses
strdrug=strcat('Pulse d=1/',num2str(pulse),' drug=',num2str((drug-1)/0.5));
title(strdrug,'FontName','Arial', 'FontSize', 28); 

%%use for different clock parameters
%strclock=strcat('with Clock A=',num2str(clockamp),' Pcd=',num2str(pcd)); 
%title(strclock,'FontName','Arial', 'FontSize', 28);

%%use for raw ppERK plots
%title('Raw ppERK','FontName','Arial', 'FontSize', 28); 

axes1=gca;
xlim(axes1,[0 400]);
ylim(axes1,[0 35]);
zlim(axes1,[-1 40]);
view(90,90);
% Set the remaining axes properties
set(gca,'BoxStyle','full','FontName','Arial','FontSize',20,'LineWidth',...
    1.5,'XGrid','on','XTick',[0 70 140 210 280 350],'YGrid','on','YTick',...
    [0 5 10 15 20 25 30 35 40]);
xticklabels({'0','30','60','90','120','150'});
caxis([0 1]);
% Create colorbar
colorbar(gca,'FontSize',16);
colorbar('Ticks',[0,0.25,0.5,0.75,1],'TickLabels',{'0','0.11','0.22','0.33','0.44'}); %use for rescaled SFC plots
%end
