function res=showMatrixAbs(cmd,varargin)
% showMatrixAbs: show absolute value of matrix entries as colored image
% function res=showMatrixAbs(cmd,varargin)
% 
% res = showMatrixAbs('init',opt)
%   Initialize figure and Color-Levels
%     opt  struct with following fields
%          fig:      (optional) figure-handle to use 
%          levels:   vector with different absolute-values where
%                    color changes
%     res  result; must be used for further calls
%
% res = showMatrix('update',res,A)
%   Update figure for new matrix A
%

switch cmd
    case 'init'
        res=initFigure(varargin{1});
    case 'update'
        res=updateFigure(varargin{:});
    otherwise
        error(['Unknown command "',cmd,'"']);
end

end

function res=initFigure(opt)
res=struct('fig',[],'ax1',[],'levels',[],'ih',[],'cb',[],'diag_lh',[]);
fig = [];
if isfield(opt,'fig'), fig=opt.fig; end
if isempty(fig), fig=figure();end
clf(fig);
levels=sort(opt.levels);
cm=jet(numel(levels)+1);
colormap(fig,cm);
ax1=axes('Parent',fig,'Box','on',...
    'Units','Normalized','Position',[0.05,0.05,0.85,0.9]);
C=uint8(0);
ih=image(C,'Parent',ax1);
set(ih,'Visible','off');
hold on;
res.diag_lh=plot(ax1,NaN,NaN,'LineWidth',2,'Color','k');hold off;
set(ax1,'XTick',[],'YTick',[],'YDir','reverse');
res.cb=colorbar(ax1,'eastoutside','AxisLocation','out',...
    'Ticks',1:numel(levels),...
    'TickLabels',cellfun(@num2str,num2cell(levels),'UniformOutput',0));
res.fig=fig;res.ax1=ax1;res.levels=levels;res.ih=ih;
end

function res=updateFigure(res,A)
levels=res.levels;
M=abs(A); C=uint8(numel(levels)+1)*ones(size(M),'uint8');
C(M<levels(1))=1;
for k=2:numel(levels)
    C( levels(k-1)<=M & M<levels(k) ) = uint8(k);
end
[m,n]=size(C);
set(res.ih,'XData',[1,n],'YData',[1,m],...
    'CData',C,'CDataMapping','direct','Visible','on');
set(res.ax1,'XLim',[0.5,n+0.5],'YLim',[0.5,m+0.5]);
axis(res.ax1,'equal');
h=kron(0:size(A,1),[1,1]);
set(res.diag_lh,'XData',h(1:end-1)+0.5,'YData',h(2:end)+0.5);
end
