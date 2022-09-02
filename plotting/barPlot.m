
function [xCoordinates,lgdObject,handleObject] = barPlot(inputData,NameValueArgs)
% barPlot Plot barplot from vector or matrix input.
%
% Created by Ryan Gorzek
%
% Dependencies: none
% 
% Input Arguments:
%
%   inputData -- vector or matrix of data for barplot. If inputData is a vector, specify
%   categorical name-value argument inputLabels to plot multiple bars. If inputData is
%   a matrix, each bar represents a column.
%
%   Name-Value Arguments (default):
%
%     inputLabels ([]) -- vector of categorical labels for vector input.
%
%     groupSize (1) -- scalar that specifies the number of bars by which to group the data (if grouping).
%
%     labelGroup (false) -- logical that specifies whether to shrink the number of x-axis labels to one per group (if grouping).
%
%     barLabels ([1:nBars]) -- cell array of strings that specify x-axis labels.
%
%     barColors ([0.7,0.7,0.7]) -- cell array of RGB vectors that specify bar colors.
%
%     barEdgeColors ([0,0,0]) -- 
%
%     barLineWidth (1) -- scalar that specifies the line width of the bar edges and whiskers.
%
%     whiskerColors ([0,0,0]) -- 
%
%     outlierSize (30) -- scalar that specifies the size of outlier points.
%
%     plotPoints (false) -- 
%
%     pointSize (10) -- 
%
%     pointColor ([0,0,0]) -- 
%
%     connectGroups (false) -- 
%
%     connectLineWidth ([0,0,0]) -- 
%
%     connectLineColor ([0,0,0]) -- 
%
%     axFontSize (13) -- scalar that specifies the axes font size.
%
%     plotLegend (false) -- logical that specifies whether or not to plot a legend with the specified parameters.
%
%     lgdLabels (none) -- cell array of strings that specify legend labels.
%
%     lgdColors (none) -- cell array of RBG vectors that specify legend colors.
%
%     lgdLocation ('northeast') -- string that specifies legend location. See MATLAB legend documentation for options.
%
%     lgdbar ('off') -- string that specifies whether to show legend bar. Options are 'on' or 'off'.
%
%     lgdFontSize (12) -- scalar that specifies legend font size.
%
%     lgdLineWidth (8) -- scalar that specifies the width of legend markers.
%
% Output Arguments:
%
%   xCoors -- vector of x-axis coordinates corresponding to each plotted bar.
%
%   lgdObject -- legend object. See MATLAB legend documentation for more information.
%
%   handleObject -- legend handle object, i.e. a graphics array of legend lines and text.
%

arguments
    
    inputData double
    NameValueArgs.inputLabels (:,1) double = []
    NameValueArgs.groupSize (1,1) {mustBeNumeric} = 1
    NameValueArgs.labelGroup = false
    NameValueArgs.barLabels (1,:) {mustBeA(NameValueArgs.barLabels,'cell'),mustBeText} = {}
    NameValueArgs.barColors (1,:) {mustBeA(NameValueArgs.barColors,'cell')} = {}
    NameValueArgs.barEdgeColors (1,:) {mustBeA(NameValueArgs.barEdgeColors,'cell')} = {}
    NameValueArgs.barLineWidth (1,1) {mustBeNumeric} = 1
    NameValueArgs.errorType (1,1) string = 'SEM'
    NameValueArgs.whiskerColors (1,:) {mustBeA(NameValueArgs.whiskerColors,'cell')} = {}
    % NameValueArgs.outlierSize (1,1) double = 30
    NameValueArgs.plotPoints = false
    NameValueArgs.pointSize (1,1) double = 10
    NameValueArgs.pointColor (1,3) double = [0,0,0]
    NameValueArgs.connectGroups = false
    NameValueArgs.connectLineWidth (1,1) double = 0.5
    NameValueArgs.connectLineColor (1,3) double = [0,0,0]
    NameValueArgs.axFontSize (1,1) double = 13
    NameValueArgs.plotLegend = true
    NameValueArgs.lgdLabels (1,:) {mustBeA(NameValueArgs.lgdLabels,'cell'),mustBeText} = {}
    NameValueArgs.lgdColors (1,:) {mustBeA(NameValueArgs.lgdColors,'cell')} = {}
    NameValueArgs.lgdLocation (1,1) string = 'northeast'
    NameValueArgs.lgdPosition (1,4) double = [0,0,0,0]
    NameValueArgs.lgdBox (1,1) string = 'off'
    NameValueArgs.lgdFontSize (1,1) double = 12
    NameValueArgs.lgdLineWidth (1,1) double = 8

end

%%%% assign temporary variables for name-value arguments

inputLabels = NameValueArgs.inputLabels;
groupSize = NameValueArgs.groupSize;
labelGroup = NameValueArgs.labelGroup;
barLabels = NameValueArgs.barLabels;
barColors = NameValueArgs.barColors;
barEdgeColors = NameValueArgs.barEdgeColors;
barLineWidth = NameValueArgs.barLineWidth;
errorType = NameValueArgs.errorType;
whiskerColors = NameValueArgs.whiskerColors;
% outlierSize = NameValueArgs.outlierSize;
plotPoints = NameValueArgs.plotPoints;
pointSize = NameValueArgs.pointSize;
pointColor = NameValueArgs.pointColor;
connectGroups = NameValueArgs.connectGroups;
connectLineWidth = NameValueArgs.connectLineWidth;
connectLineColor = NameValueArgs.connectLineColor;
axFontSize = NameValueArgs.axFontSize;
plotLegend = NameValueArgs.plotLegend;
lgdLabels = NameValueArgs.lgdLabels;
lgdColors = NameValueArgs.lgdColors;
lgdLocation = NameValueArgs.lgdLocation;
lgdPosition = NameValueArgs.lgdPosition;
lgdBox = NameValueArgs.lgdBox;
lgdFontSize = NameValueArgs.lgdFontSize;
lgdLineWidth = NameValueArgs.lgdLineWidth;

%%%% check for vector or matrix input data

if isempty(inputLabels) % if no inputLabels specified, reshape matrix input into vector and produce labels
    
    inputLabels = repmat(1:size(inputData,2),[size(inputData,1),1]); % create labels for inputData matrix
    
    inputLabels = reshape(inputLabels,[],1); % reshape inputLabels into a column vector
    inputData = reshape(inputData,[],1);    % reshape inputData into a column vector
    
elseif ~isempty(inputLabels) && size(inputData,2) ~= 1 % if inputLabels are specified but inputData is not in column vector, reshape
    
    inputData = reshape(inputData,[],1);

end

%%%% calculate the number of bars and groups

nBars = numel(unique(inputLabels)); % get number of bars to plot based on inputLabels
nGroups = nBars/groupSize; % get number of bar groups based on nBars and groupSize

uniqueLabels = sort(unique(inputLabels),'ascend')'; % get unique labels from inputLabels

%%%% throw error if nGroups is not an integer

if rem(nGroups,1) ~= 0, error('Number of input categories is not divisible by number of groups.'); end

%%%% set default bar labels (numbered) and bar colors (gray) if not specified

if isempty(barLabels) && groupSize == 1, barLabels = cellstr(string(1:nBars)); elseif isempty(barLabels) && groupSize > 1 && labelGroup == false, barLabels = cellstr(string(1:groupSize)); elseif isempty(barLabels) && labelGroup == true, barLabels = cellstr(string(1:nGroups)); end
if isempty(barColors), barColors = repmat({[0.7,0.7,0.7]},[1,numel(unique(inputLabels))]); elseif numel(barColors) == groupSize, barColors = repmat(barColors,[1,nGroups]); end
if isempty(barEdgeColors), barEdgeColors = repmat({[0,0,0]},[1,numel(unique(inputLabels))]); elseif numel(barEdgeColors) == groupSize, barEdgeColors = repmat(barEdgeColors,[1,nGroups]); end
if isempty(whiskerColors), whiskerColors = repmat({[0,0,0]},[1,numel(unique(inputLabels))]); elseif numel(whiskerColors) == groupSize, whiskerColors = repmat(whiskerColors,[1,nGroups]); end

%%%% check barLabels against nGroups to warn user about specifying labelGroup

if labelGroup == true && numel(barLabels) > nGroups
    error('Too many barLabels for number of groups, did you mean to specify labelGroup = true?');
elseif labelGroup == true && numel(barLabels) < nGroups
    error('Too few barLabels for number of groups, did you mean to specify labelGroup = true?');
end

%%%% check whether labels match nGroups if labelGroup = false

if groupSize ~= 1 && ...
   numel(barLabels) == nGroups && ...
   numel(barColors) == groupSize && ...
   labelGroup == false
    
    error('Insufficient number of barLabels, did you mean to specify labelGroup = true?');

end

%%%% get x-coordinates for plotting

if nGroups == 1
    
    xCoordinates = [0.70:0.65:0.7+(0.65*nBars)-0.65];
   
else

    xCoordinates = [0.70:0.55:0.70+(0.55*(nBars/nGroups)-0.55)]; initCoors = [0.70:0.55:0.70+(0.55*(nBars/nGroups)-0.55)];

    for grp = 2:nGroups, xCoordinates = horzcat(xCoordinates,initCoors+(xCoordinates(end)+0.4)); end

end

%%%% intialize matrices for storing max/min to set axes

maxMat = zeros(nBars,3); % uQuar, uWhisk, outliersU
minMat = zeros(nBars,3); % lQuar, lWhisk, outliersL

%%%% plot

for bar = uniqueLabels

    barNum = find(uniqueLabels == bar);
    
    clear currMean

    %%%% get data for current bar
    
    currData = inputData(inputLabels == bar,1);
    
    %%%% plot bar if there is at least one data point
    
    if nnz(~isnan(currData)) > 0
        
        % get mean of categories
        currMean = mean(currData,1,'omitnan');

        if ~isempty(errorType) && strcmp(errorType,'STD'), errorBounds = std(currData,0,1,'omitnan'); end

        if ~isempty(errorType) && strcmp(errorType,'SEM'), errorBounds = std(currData,0,1,'omitnan')./sqrt(sum(~isnan(currData))); end
        
        hold on;

        % plot bar with bounds at 0 and mean
        if currMean >= 0, rectangle('Position',[xCoordinates(barNum)-0.25 0 0.5 currMean],'FaceColor',barColors{barNum},'EdgeColor',barEdgeColors{barNum},'LineWidth',barLineWidth,'Tag','bar');

        else, rectangle('Position',[xCoordinates(barNum)-0.25 currMean 0.5 -1*currMean],'FaceColor',barColors{barNum},'EdgeColor',barEdgeColors{barNum},'LineWidth',barLineWidth,'Tag','bar');

        end

        % plot lower whisker
        line([xCoordinates(barNum) xCoordinates(barNum)],[currMean,currMean - errorBounds],'Color',whiskerColors{barNum},'LineWidth',barLineWidth);
        % plot upper whisker
        line([xCoordinates(barNum) xCoordinates(barNum)],[currMean,currMean + errorBounds],'Color',whiskerColors{barNum},'LineWidth',barLineWidth,'Tag','uwhisk');
        % plot lower whisker bar
        line([xCoordinates(barNum)-0.1 xCoordinates(barNum)+0.1],[currMean - errorBounds currMean - errorBounds],'Color',whiskerColors{barNum},'LineWidth',barLineWidth);
        % plot upper whisker bar
        line([xCoordinates(barNum)-0.1 xCoordinates(barNum)+0.1],[currMean + errorBounds currMean + errorBounds],'Color',whiskerColors{barNum},'LineWidth',barLineWidth);

        %%%% plot outliers (?)

        %%%% store min/max from each category for setting axis limits
        maxMat(barNum,:) = [0,currMean - errorBounds,currMean + errorBounds];
        minMat(barNum,:) = [0,currMean - errorBounds,currMean + errorBounds];
        
    elseif nnz(~isnan(currData)) == 0
        
        %%%% store min/max from each category for setting limits
        maxMat(barNum,:) = [repmat(max(currData),[1,3])];
        minMat(barNum,:) = [repmat(min(currData),[1,3])];
        
    end
    
end

%%%% plot individual data points and/or connect groups

if (isnumeric(plotPoints) || plotPoints == true) && nGroups == nBars

    currData = reshape(inputData,[],nBars);

    if isnumeric(plotPoints), pointIdx = plotPoints; else, pointIdx = 1:nBars; end

    for bar = pointIdx
        
        scatter(xCoordinates(bar),currData(:,bar),pointSize,'MarkerFaceColor',pointColor,'MarkerEdgeColor',pointColor);

    end

elseif (isnumeric(plotPoints) || plotPoints == true) && nGroups ~= nBars

    currData = reshape(inputData,[],nBars);

    if isnumeric(plotPoints), pointIdx = plotPoints; else, pointIdx = 1:nBars; end

    barIdx = reshape(1:nBars,[],nGroups)'; barIdx = barIdx(:,pointIdx);
    
    for group = 1:nGroups

        for bar = 1:size(barIdx,2), scatter(xCoordinates(barIdx(group,bar)),currData(:,barIdx(group,bar)),pointSize,'MarkerFaceColor',pointColor,'MarkerEdgeColor',pointColor); end

    end

end

if (isnumeric(connectGroups) || connectGroups == true) && nGroups == nBars

    currData = reshape(inputData,[],nBars);

    if isnumeric(connectGroups), connectIdx = connectGroups; else, connectIdx = 1:nBars; end

    for sample = 1:size(currData,1), plot(xCoordinates(connectIdx),currData(sample,connectIdx),'LineWidth',connectLineWidth,'Color',connectLineColor); end

elseif (isnumeric(connectGroups) || connectGroups == true) && nGroups ~= nBars

    currData = reshape(inputData,[],nBars);

    if isnumeric(connectGroups), connectIdx = connectGroups; else, connectIdx = 1:nBars; end
    
    for group = 1:nGroups

        barIdx = reshape(1:nBars,[],nGroups)'; barIdx = barIdx(:,connectIdx);
        
        for sample = 1:size(currData,1), plot(xCoordinates(barIdx(group,:)),currData(sample,barIdx(group,:)),'LineWidth',connectLineWidth,'Color',connectLineColor); end

    end

end

%%%% set x-axis limits

xlim([0.2 xCoordinates(end)+0.5]);

%%%% set x-axis tick labels

xLabPos = [];

if labelGroup == true

    for group = 1:(nBars/nGroups):nBars, xLabPos = horzcat(xLabPos,median(xCoordinates(group:group+(nBars/nGroups)-1))); end

    set(gca,'xtick',xLabPos); set(gca,'xticklabel',barLabels);

else
        
    set(gca,'xtick',xCoordinates); set(gca,'xticklabel',barLabels);

end

%%%% set figure & axes appearance

set(gcf,'color','w');
set(gca,'box','off','XColor','k','YColor','k','TickDir','out','TickLength',[0.01,0.01],'FontSize',axFontSize,'LineWidth',1);

%%%% set y-axis limits

if ~all(isnan(maxMat),'all')

    yUpper = max(maxMat,[],'all'); yLower = min(minMat,[],'all');
    
    yExt = (yUpper-yLower)*0.2;
    
    if yLower < 0, ylim([yLower-yExt yUpper+yExt]); else, ylim([0,yUpper+yExt]); end

end

%%%% add legend if specified

if plotLegend == true && ~isempty(lgdLabels) && ~isempty(lgdColors)

    [lgdObject,handleObject,~,~] = legend(lgdLabels);
    
    h1 = findobj(handleObject,'Type','line'); set(h1,'LineWidth',lgdLineWidth);

    lineLoc = 1; for col = 1:numel(lgdColors), set(h1(lineLoc),'Color',lgdColors{col}); lineLoc = lineLoc + 2; end

    if strcmp(lgdBox,'off'), legend boxoff; end

    lgdObject.Location = lgdLocation;
    if any(lgdPosition), lgdObject.Position = lgdPosition; end
    lgdObject.String = strcat('\fontsize{',num2str(lgdFontSize),'}',lgdObject.String);
    lgdObject.FontSize = lgdFontSize;
    
else
    
    lgdObject = []; handleObject = [];

end

end
