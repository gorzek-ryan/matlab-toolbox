
function setStyle(NameValueArgs)
% setStyle Set axes/figure style, labels, limits, and position.
%
% Author: Ryan Gorzek
%
% Dependencies: none
%
% Input Arguments:
%
%   Name-Value Arguments (default):
%
%     title (none) -- string that specifies plot title.
%
%     xlabel (none) -- string that specifies x-axis label.
%
%     xticks (automatic) -- vector that specifies x-axis tick marks.
%
%     xticklabels (automatic) -- cell array of strings that specify x-axis tick labels.
%
%     xtickangle (automatic) -- scalar that specifies x-axis tick label angle.
%
%     xlim (automatic) -- vector that specifies x-axis limits.
%
%     ylabel (none) -- string that specifies y-axis label.
%
%     yticks (automatic) -- vector that specifies y-axis tick marks.
%
%     yticklabels (automatic) -- cell array of strings that specify y-axis tick labels.
%
%     ytickangle (automatic) -- scalar that specifies y-axis tick label angle.
%
%     ylim (automatic) -- vector that specifies y-axis limits.
%
%     fontSize (13) -- scalar that specifies font size of axes text.
%
%     lineWidth (1) -- scalar that specifies axes line width.
%
%     figPosition (automatic) -- vector that specifies figure position.
%

arguments
    
    NameValueArgs.title (1,1) string = ''
    NameValueArgs.xlabel (1,1) string = ''
    NameValueArgs.xticks (1,:) double = []
    NameValueArgs.xticklabels = []
    NameValueArgs.xtickangle = []
    NameValueArgs.xlim = []
    NameValueArgs.ylabel (1,1) string = ''
    NameValueArgs.yticks (1,:) double = []
    NameValueArgs.yticklabels = []
    NameValueArgs.ytickangle = []
    NameValueArgs.ylim = []
    NameValueArgs.fontSize (1,1) double = 13
    NameValueArgs.lineWidth (1,1) double = 1
    NameValueArgs.figPosition = []

end

if ~isempty(NameValueArgs.title), title(NameValueArgs.title); end
if ~isempty(NameValueArgs.xlabel), xlabel(NameValueArgs.xlabel); end
if isnan(NameValueArgs.xticks), xticks([]); elseif ~isempty(NameValueArgs.xticks), xticks(NameValueArgs.xticks); end
if ~isempty(NameValueArgs.xticklabels), xticklabels(NameValueArgs.xticklabels); end
if ~isempty(NameValueArgs.xtickangle), xtickangle(NameValueArgs.xtickangle); end
if ~isempty(NameValueArgs.xlim), xlim(NameValueArgs.xlim); end
if ~isempty(NameValueArgs.ylabel), ylabel(NameValueArgs.ylabel); end
if isnan(NameValueArgs.yticks), yticks([]); elseif ~isempty(NameValueArgs.yticks), yticks(NameValueArgs.yticks); end
if isnan(NameValueArgs.yticklabels), yticklabels([]); elseif ~isempty(NameValueArgs.yticklabels), yticklabels(NameValueArgs.yticklabels); end
if ~isempty(NameValueArgs.ytickangle), ytickangle(NameValueArgs.ytickangle); end
if ~isempty(NameValueArgs.ylim), ylim(NameValueArgs.ylim); end

fontSize = NameValueArgs.fontSize;
lineWidth = NameValueArgs.lineWidth;

set(gcf,'color','w');
set(gca,'box','off','XColor','k','YColor','k','TickDir','out','TickLength',[0.01,0.01],'FontSize',fontSize,'LineWidth',lineWidth);

if ~isempty(NameValueArgs.figPosition), set(gcf,'Position',NameValueArgs.figPosition); end

end
