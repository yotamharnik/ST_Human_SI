function plot_patch(x,mn,se,color,PATCH_ONLY,LineW)

% Plots a patch with mean+- one standard error and the given color
% PTCH_ONLY=0 - draw both, 1 - draw only patch, 2 - draw only line 
% x=0:20;
% y1=smooth(mn1(21:end));
% y2=smooth(se1(21:end));

mn(isnan(mn))=0;
se(isnan(se))=0;
if nargin<5,
    PATCH_ONLY=0;
end
if nargin<6,
    LineW=2;
end
x=x(:);mn=mn(:);se=se(:);
hold on;
if ischar(color)
    if contains(color,'#')
        color2=hex2rgb(color);
    else
        color2=color;
    end
else
    color2=color;
end

if PATCH_ONLY==0,
    plot(x,mn,'LineWidth',LineW,'color',color);
    patch([x ;x(end:-1:1)],[mn+se; mn(end:-1:1)-se(end:-1:1)],color2,'EdgeColor','none','FaceAlpha',0.2);
elseif PATCH_ONLY==1,
    patch([x ;x(end:-1:1)],[mn+se; mn(end:-1:1)-se(end:-1:1)],color2,'EdgeColor','none','FaceAlpha',0.2);
else
    plot(x,mn,'LineWidth',LineW,'color',color);
end
%axis square;
%axis([xlim 0 1.1]);
set(gca,'FontSize',16);
set(gca,'xtick',1:length(x));
%line([3 3],ylim,'color','k','linestyle','--','linewidth',2);
%line([6 6],ylim,'color','k','linestyle','--','linewidth',2);
axis tight;