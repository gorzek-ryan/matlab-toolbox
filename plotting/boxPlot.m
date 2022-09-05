
function [xCoordinates,lgdObject] = boxPlot(inputData,NameValueArgs)
% boxPlot Plot boxplot from vector or matrix input.
%
% Created by Ryan Gorzek
%
% Dependencies: none
% 
% Input Arguments:
%
%   inputData -- vector or matrix of data for boxplot. If inputData is a matrix, 
%                each column will be plotted as a box.If inputData is a vector, specify
%                categorical name-value argument inputLabels to plot multiple boxes. 
%
%   Name-Value Arguments (default):
%
%     inputLabels ([]) -- vector of categorical labels for vector input.
%
%     groupSize (1) -- scalar that specifies the number of boxes by which to group the data (if grouping).
%
%     labelGroup (false) -- logical that specifies whether to shrink the number of x-axis labels to one per group (if grouping).
%
%     boxLabels ([1:nBoxes]) -- cell array of strings that specify x-axis labels.
%
%     boxColors ([0.7,0.7,0.7]) -- cell array of RGB vectors that specify box colors.
%
%     boxEdgeColors ([0,0,0]) -- 
%
%     boxLineWidth (1) -- scalar that specifies the line width of the box edges and whiskers.
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
%     plotLegend (false) -- logical that specifies whether or not to plot a legend with the specified parameters.
%
%     lgdLabels (none) -- cell array of strings that specify legend labels.
%
%     lgdColors (none) -- cell array of RBG vectors that specify legend colors.
%
%     lgdColumns (1) -- 
%
%     lgdOrientation ('horizontal') -- 
%
%     lgdBox ('off') -- string that specifies whether to show legend box. Options are 'on' or 'off'.
%
%     lgdFontSize (12) -- scalar that specifies legend font size.
%
%     lgdLineHeight (1) -- scalar that specifies the height of legend markers.
%
%     lgdLineWidth (1) -- scalar that specifies the width of legend markers.
%
%     lgdLocation ('northeast') -- string that specifies legend location. See MATLAB legend documentation for options.
%
%     lgdPosition ([]) -- 
%
% Output Arguments:
%
%   xCoors -- vector of x-axis coordinates corresponding to each plotted box.
%
%   lgdObject -- legend object. See MATLAB legend documentation for more information.
%

arguments
    
    inputData (:,:) {mustBeNumeric} % 
    
    NameValueArgs.inputLabels {mustBeVector,mustBeNumeric} = reshape(repmat(1:size(inputData,2),[size(inputData,1),1]),[],1) % 

    NameValueArgs.groupSize (1,1) {mustBeNumeric} = 1 % 
    NameValueArgs.labelGroups (1,1) logical = false % 

    NameValueArgs.boxLabels {mustBeA(NameValueArgs.boxLabels,'cell')} = {} %
    NameValueArgs.boxColors {mustBeA(NameValueArgs.boxColors,'cell')} = {} % 
    NameValueArgs.boxAlpha (1,1) {mustBeInRange(NameValueArgs.boxAlpha,0,1)} = 1
    NameValueArgs.boxEdgeColors {mustBeA(NameValueArgs.boxEdgeColors,'cell')} = {} %
    NameValueArgs.boxEdgeWidth (1,1) {mustBeNumeric} = 1 %
    NameValueArgs.boxEdgeStyle (1,1) string = '-' %
    NameValueArgs.boxEdgeAlpha (1,1) {mustBeInRange(NameValueArgs.boxEdgeAlpha,0,1)} = 1
    NameValueArgs.boxSpacing (1,1) {mustBeNumeric} = 1 %
    NameValueArgs.boxOrientation {mustBeMember(NameValueArgs.boxOrientation,['vertical','horizontal'])} = 'vertical' %%%% incorporate this
    NameValueArgs.boxCurvature (1,2) double = [0,0] % 

    NameValueArgs.medianColors (1,:) {mustBeA(NameValueArgs.medianColors,'cell')} = {} % 
    NameValueArgs.medianWidth (1,1) {mustBeNumeric} = 2 % 
    NameValueArgs.medianStyle (1,1) string = '-'  % 
    NameValueArgs.medianAlpha (1,1) {mustBeInRange(NameValueArgs.medianAlpha,0,1)} = 1

    NameValueArgs.whiskerColors (1,:) {mustBeA(NameValueArgs.whiskerColors,'cell')} = {} % 
    NameValueArgs.whiskerWidth (1,1) {mustBeNumeric} = 1 % 
    NameValueArgs.whiskerStyle (1,1) string = '-' % 
    NameValueArgs.whiskerAlpha (1,1) {mustBeInRange(NameValueArgs.whiskerAlpha,0,1)} = 1

    NameValueArgs.outlierColors (1,:) {mustBeA(NameValueArgs.outlierColors,'cell')} = {} % 
    NameValueArgs.outlierSize (1,1) {mustBeNumeric} = 30 % 
    NameValueArgs.outlierStyle (1,1) string = 'o' % 
    NameValueArgs.outlierAlpha (1,1) {mustBeInRange(NameValueArgs.outlierAlpha,0,1)} = 1
    NameValueArgs.outlierJitter {mustBeMember(NameValueArgs.outlierJitter,['none','density','rand','randn'])} = 'none'

    NameValueArgs.plotPoints {mustBeNumericOrLogical} = false
    NameValueArgs.pointColors {mustBeA(NameValueArgs.pointColors,'cell')} = {}
    NameValueArgs.pointSize (1,1) {mustBeNumeric} = 20 % 
    NameValueArgs.pointStyle (1,1) string = '.' % 
    NameValueArgs.pointAlpha (1,1) {mustBeInRange(NameValueArgs.pointAlpha,0,1)} = 1
    NameValueArgs.pointJitter {mustBeMember(NameValueArgs.pointJitter,['none','density','rand','randn'])} = 'none'

    NameValueArgs.plotLines {mustBeNumericOrLogical} = false
    NameValueArgs.lineColors {mustBeA(NameValueArgs.lineColors,'cell')} = {}
    NameValueArgs.lineWidth (1,1) {mustBeNumeric} = 0.5 % 
    NameValueArgs.lineStyle (1,1) string = '-' % 
    NameValueArgs.lineAlpha (1,1) {mustBeInRange(NameValueArgs.lineAlpha,0,1)} = 1

    NameValueArgs.plotLegend (1,1) logical = false % 
    NameValueArgs.lgdLabels (1,:) {mustBeA(NameValueArgs.lgdLabels,'cell')} = {} % 
    NameValueArgs.lgdColors (1,:) {mustBeA(NameValueArgs.lgdColors,'cell')} = {} % 
    NameValueArgs.lgdColumns (1,1) {mustBeNumeric} = 1 % 
    NameValueArgs.lgdOrientation {mustBeMember(NameValueArgs.lgdOrientation,['vertical','horizontal'])} = 'horizontal'
    NameValueArgs.lgdBox {mustBeMember(NameValueArgs.lgdBox,['on','off'])} = 'off'
    % background color
    % box color
    NameValueArgs.lgdFontSize (1,1) {mustBeNumeric} = 12 % 
    NameValueArgs.lgdFontWeight {mustBeMember(NameValueArgs.lgdFontWeight,['normal','bold'])} = 'normal' %%%%
    % font color
    % font angle?
    NameValueArgs.lgdLineHeight (1,1) {mustBeNumeric} = 1 % 
    NameValueArgs.lgdLineWidth (1,1) {mustBeNumeric} = 1 % 
    NameValueArgs.lgdLocation (1,1) string = 'northeast' % 
    NameValueArgs.lgdPosition (1,4) {mustBeInRange(NameValueArgs.lgdPosition,0,1)} = [0,0,0,0] % 

end

% Get number of boxes, groups, samples (including NaN), and missing values
nBoxes = numel(unique(NameValueArgs.inputLabels));
nGroups = nBoxes/NameValueArgs.groupSize;
nSamples = size(inputData,1);
nMissing = nnz(all(isnan(inputData),2));

% Throw error if number of groups is not an integer
if rem(nGroups,1) ~= 0
    error('Number of input boxes is not divisible by number of groups.'); 
end

% Reshape input data into a column vector %%%% and get unique labels
inputData = reshape(inputData,[],1);
uniqueLabels = sort(unique(NameValueArgs.inputLabels),'ascend')'; %%%%

% Set default box labels (numbered) if none are specified
if isempty(NameValueArgs.boxLabels) && ...
   NameValueArgs.groupSize == 1

    NameValueArgs.boxLabels = cellstr(string(1:nBoxes));

elseif isempty(NameValueArgs.boxLabels) && ...
       NameValueArgs.groupSize > 1 && ...
       NameValueArgs.labelGroups == false

    NameValueArgs.boxLabels = cellstr(string(1:NameValueArgs.groupSize));

elseif isempty(NameValueArgs.boxLabels) && ...
       NameValueArgs.labelGroups == true
    
    NameValueArgs.boxLabels = cellstr(string(1:nGroups));

end

% Set default box colors (gray) if none are specified, or replicate
% single or group-level color specification
if isempty(NameValueArgs.boxColors)
    NameValueArgs.boxColors = repmat({[0.7,0.7,0.7]},[1,nBoxes]);
elseif numel(NameValueArgs.boxColors) == 1
    NameValueArgs.boxColors = repmat(NameValueArgs.boxColors,[1,nBoxes]);
elseif numel(NameValueArgs.boxColors) == groupSize
    NameValueArgs.boxColors = repmat(NameValueArgs.boxColors,[1,nGroups]);
end

% Set default box edge colors (black) if none are specified, or replicate
% single or group-level color specification
if isempty(NameValueArgs.boxEdgeColors)
    NameValueArgs.boxEdgeColors = repmat({[0.0,0.0,0.0]},[1,nBoxes]);
elseif numel(NameValueArgs.boxEdgeColors) == 1
    NameValueArgs.boxEdgeColors = repmat(NameValueArgs.boxEdgeColors,[1,nBoxes]);
elseif numel(NameValueArgs.boxEdgeColors) == groupSize
    NameValueArgs.boxEdgeColors = repmat(NameValueArgs.boxEdgeColors,[1,nGroups]);
end

% Set default whisker colors (black) if none are specified, or replicate
% single or group-level color specification
if isempty(NameValueArgs.whiskerColors)
    NameValueArgs.whiskerColors = repmat({[0.0,0.0,0.0]},[1,nBoxes]);
elseif numel(NameValueArgs.whiskerColors) == 1
    NameValueArgs.whiskerColors = repmat(NameValueArgs.whiskerColors,[1,nBoxes]);
elseif numel(NameValueArgs.whiskerColors) == groupSize
    NameValueArgs.whiskerColors = repmat(NameValueArgs.whiskerColors,[1,nGroups]);
end

% Set default median colors (black) if none are specified, or replicate
% single or group-level color specification
if isempty(NameValueArgs.medianColors)
    NameValueArgs.medianColors = repmat({[0.0,0.0,0.0]},[1,nBoxes]);
elseif numel(NameValueArgs.medianColors) == 1
    NameValueArgs.medianColors = repmat(NameValueArgs.medianColors,[1,nBoxes]);
elseif numel(NameValueArgs.medianColors) == groupSize
    NameValueArgs.medianColors = repmat(NameValueArgs.medianColors,[1,nGroups]);
end

% Set default outlier colors (matched to box colors) if none are specified,
% or replicate single or group-level color specification
if isempty(NameValueArgs.outlierColors)
    NameValueArgs.outlierColors = NameValueArgs.boxColors;
elseif numel(NameValueArgs.outlierColors) == 1
    NameValueArgs.outlierColors = repmat(NameValueArgs.outlierColors,[1,nBoxes]);
elseif numel(NameValueArgs.outlierColors) == groupSize
    NameValueArgs.outlierColors = repmat(NameValueArgs.outlierColors,[1,nGroups]);
end

% Set default point colors (black) if not specified ...
%
if isempty(NameValueArgs.pointColors)
    NameValueArgs.pointColors = repmat({[0.0,0.0,0.0]},[nSamples,nBoxes]);
elseif numel(NameValueArgs.pointColors) == 1
    NameValueArgs.pointColors = repmat(NameValueArgs.pointColors,[nSamples,nBoxes]);
elseif isvector(NameValueArgs.pointColors) && ...
       any(size(NameValueArgs.pointColors) == groupSize)
    NameValueArgs.pointColors = repmat(NameValueArgs.pointColors,[nSamples,nGroups]);
elseif isvector(NameValueArgs.pointColors) && ...
       any(size(NameValueArgs.pointColors) == nBoxes)
    NameValueArgs.pointColors = repmat(NameValueArgs.pointColors,[nSamples,1]);
elseif all(ismember(size(NameValueArgs.pointColors),[nSamples - nMissing,groupSize])) %%%%
    NameValueArgs.pointColors = repmat(NameValueArgs.pointColors,[1,nGroups]);
% need condition for pointColors matching nSamples
end

% Check default color stream specification for point colors
point_defaultIdx = cell2mat(cellfun(@(x) strcmp(x,'default'),NameValueArgs.pointColors,'UniformOutput',false));
if any(point_defaultIdx,'all')
    defaultColors = repmat(num2cell(colororder,2),[ceil(size(NameValueArgs.pointColors,1)/7),size(NameValueArgs.pointColors,2)]);
    defaultColors = defaultColors(1:size(NameValueArgs.pointColors,1),:);
    NameValueArgs.pointColors(point_defaultIdx) = defaultColors(point_defaultIdx);
end

% Set default line colors (black) if not specified ...
%
nLines = nGroups*(groupSize - 1);
if isempty(NameValueArgs.lineColors)
    NameValueArgs.lineColors = repmat({[0.0,0.0,0.0]},[nSamples,nLines]);
elseif numel(NameValueArgs.lineColors) == 1
    NameValueArgs.lineColors = repmat(NameValueArgs.lineColors,[nSamples,nLines]);
elseif isvector(NameValueArgs.lineColors) && ...
       any(size(NameValueArgs.lineColors) == nGroups)
    NameValueArgs.lineColors = NameValueArgs.lineColors(kron(1:nGroups,ones(nSamples,groupSize-1)));
elseif isvector(NameValueArgs.lineColors) && ...
       any(size(NameValueArgs.lineColors) == nLines)
    NameValueArgs.lineColors = repmat(NameValueArgs.lineColors,[nSamples,1]);
elseif all(ismember(size(NameValueArgs.lineColors),[nSamples - nMissing,nGroups])) %%%%
    NameValueArgs.lineColors = NameValueArgs.lineColors(:,kron(1:nGroups,ones(1,groupSize-1)));
% need condition for lineColors matching nSamples
end

% Check default color stream specification for line colors
line_defaultIdx = cell2mat(cellfun(@(x) strcmp(x,'default'),NameValueArgs.lineColors,'UniformOutput',false));
if any(line_defaultIdx,'all')
    defaultColors = repmat(num2cell(colororder,2),[ceil(size(NameValueArgs.lineColors,1)/7),size(NameValueArgs.lineColors,2)]); 
    defaultColors = defaultColors(1:size(NameValueArgs.lineColors,1),:);
    NameValueArgs.lineColors(line_defaultIdx) = defaultColors(line_defaultIdx);
end

% Check boxLabels against nGroups to warn user about specifying labelGroups
if NameValueArgs.labelGroups == true && ...
   numel(NameValueArgs.boxLabels) > nGroups
    error('Number of box labels exceeds number of groups, did you mean to specify labelGroups = false?');
elseif NameValueArgs.labelGroups == true && ...
       numel(NameValueArgs.boxLabels) < nGroups
    error('Insufficient number of box labels for number of groups.');
elseif NameValueArgs.labelGroups == false && ...
       NameValueArgs.groupSize ~= 1 && ...
       numel(NameValueArgs.boxLabels) == nGroups && ...
       numel(NameValueArgs.boxColors) == NameValueArgs.groupSize %%%% why?
    error('Insufficient number of box labels for number of groups, did you mean to specify labelGroups = true?');
end

% Generate x-coordinates for plotting boxes
if nGroups == 1
    xCoordinates = (0.70 : 0.65 : 0.70 + (0.65*nBoxes) - 0.65).*boxSpacing;
else
    xCoordinates = 0.70 : 0.55 : 0.70 + (0.55*(nBoxes/nGroups) - 0.55); 
    initCoors = 0.70 : 0.55 : 0.70 + (0.55*(nBoxes/nGroups) - 0.55);
    for grp = 2:nGroups
        xCoordinates = horzcat(xCoordinates,initCoors + (xCoordinates(end) + 0.4)); %%%%
    end
    xCoordinates = xCoordinates.*boxSpacing;
end

%%%% stylize from here

% Generate matrix of jitter if specified
if strcmp(outlierJitter,'density') || strcmp(pointJitter,'density') %%%% add this

    kernelDensity = ksdensity(currData);

elseif strcmp(outlierJitter,'rand') || strcmp(pointJitter,'rand')

    jitterMat = rand(size(reshape(inputData,[],nBoxes))).*(boxWidth*0.75); jitterMat = jitterMat - mean(jitterMat,1);

elseif strcmp(outlierJitter,'randn') || strcmp(pointJitter,'randn')

    randnMat = randn(size(reshape(inputData,[],nBoxes)));

    jitterMat = (randnMat./max(randnMat,[],1)).*(boxWidth*0.75); jitterMat = jitterMat - mean(jitterMat,1);

elseif strcmp(outlierJitter,'none') && strcmp(pointJitter,'none')

    jitterMat = zeros(size(reshape(inputData,[],nBoxes)));

end

%%%% intialize matrices for storing max/min to set axes

maxMat = zeros(nBoxes,3); % upperQuantile, upperWhisker, upperOutliers
minMat = zeros(nBoxes,3); % lowerQuantile, lowerWhisker, lowerOutliers

%%%% plot

for box = uniqueLabels

    boxNum = find(uniqueLabels == box);
    
    clear boxMedian lowerQuantile upperQuantile lowerWhisker upperWhisker lowerOutliers upperOutliers

    %%%% get data for current box
    
    currData = inputData(inputLabels == box,1);
    
    %%%% plot box if there are at least 4 data points
    
    if nnz(~isnan(currData)) > 4

        boxMedian = median(currData,1,'omitnan'); % get median
        
        upperQuantile = quantile(currData,0.75); % get 75 percentile
        
        lowerQuantile = quantile(currData,0.25); % get 25 percentile

        % get upper whisker value

        maxWhisker = upperQuantile + 1.5*(upperQuantile - lowerQuantile); 
        
        upperWhisker = max(currData(currData < maxWhisker & currData >= upperQuantile));
        
        if isempty(upperWhisker), upperWhisker = maxWhisker; end

        % get lower whisker value

        minWhisker = lowerQuantile - 1.5*(upperQuantile - lowerQuantile);

        lowerWhisker = min(currData(currData > minWhisker & currData <= lowerQuantile));

        if isempty(lowerWhisker), lowerWhisker = minWhisker; end
        
        hold on;

        % plot box with bounds at quartiles
        rectangle('Position',[xCoordinates(boxNum) - 0.25,lowerQuantile,0.5,upperQuantile - lowerQuantile],...
                  'FaceColor',[boxColors{boxNum},boxAlpha],...
                  'EdgeColor',[boxEdgeColors{boxNum},boxEdgeAlpha],...
                  'LineWidth',boxEdgeWidth,...
                  'LineStyle',boxEdgeStyle,...
                  'Curvature',boxCurvature,...
                  'Tag','Box');
        
        % plot median line
        line([xCoordinates(boxNum) - 0.25,xCoordinates(boxNum) + 0.25],...
             [boxMedian,boxMedian],...
             'Color',[medianColors{boxNum},medianAlpha],...
             'LineWidth',medianWidth,...
             'LineStyle',medianStyle,...
             'Tag','Median');

        % plot lower whisker
        line([xCoordinates(boxNum),xCoordinates(boxNum)],...
             [lowerQuantile,lowerWhisker],...
             'Color',[whiskerColors{boxNum},whiskerAlpha],...
             'LineWidth',whiskerWidth,...
             'LineStyle',whiskerStyle,...
             'Tag','Lower Whisker');

        % plot upper whisker
        line([xCoordinates(boxNum),xCoordinates(boxNum)],...
             [upperQuantile,upperWhisker],...
             'Color',[whiskerColors{boxNum},whiskerAlpha],...
             'LineWidth',whiskerWidth,...
             'LineStyle',whiskerStyle,...
             'Tag','Upper Whisker');

        % plot lower whisker bar
        line([xCoordinates(boxNum) - 0.1,xCoordinates(boxNum) + 0.1],...
             [lowerWhisker,lowerWhisker],...
             'Color',[whiskerColors{boxNum},whiskerAlpha],...
             'LineWidth',whiskerWidth,...
             'LineStyle',whiskerStyle,...
             'Tag','Lower Whisker Bar');

        % plot upper whisker bar
        line([xCoordinates(boxNum) - 0.1,xCoordinates(boxNum) + 0.1],...
             [upperWhisker,upperWhisker],...
             'Color',[whiskerColors{boxNum},whiskerAlpha],...
             'LineWidth',whiskerWidth,...
             'LineStyle',whiskerStyle,...
             'Tag','Upper Whisker Bar');

        %%%% plot outliers
        
        % upper

        if any(currData > upperWhisker)
            
            % scatter outliers
            scatter(xCoordinates(boxNum) - jitterMat(currData > upperWhisker,boxNum),...
                    currData(currData > upperWhisker),...
                    outlierSize,...
                    'MarkerFaceColor',outlierColors{boxNum},...
                    'MarkerEdgeColor',outlierColors{boxNum},...
                    'Marker',outlierStyle,...
                    'MarkerFaceAlpha',outlierAlpha,...
                    'MarkerEdgeAlpha',outlierAlpha,...
                    'Tag','Outlier');

            upperOutliers = currData(currData > upperWhisker);
            
        else
            
            upperOutliers = nan;
            
        end
        
        % lower

        if any(currData < lowerWhisker)
            
            % scatter outliers
            scatter(xCoordinates(boxNum) - jitterMat(currData < lowerWhisker,boxNum),...
                    currData(currData < lowerWhisker),...
                    outlierSize,...
                    'MarkerFaceColor',outlierColors{boxNum},...
                    'MarkerEdgeColor',outlierColors{boxNum},...
                    'Marker',outlierStyle,...
                    'MarkerFaceAlpha',outlierAlpha,...
                    'MarkerEdgeAlpha',outlierAlpha,...
                    'Tag','Outlier');

            lowerOutliers = currData(currData < lowerWhisker);
            
        else
            
            lowerOutliers = nan;
            
        end

        %%%% store min/max from each box for setting axis limits

        maxMat(boxNum,:) = [upperQuantile,upperWhisker,max(upperOutliers)];
        minMat(boxNum,:) = [lowerQuantile,lowerWhisker,min(lowerOutliers)];
        
    elseif nnz(~isnan(currData)) == 0
        
        %%%% store min/max from each box for setting axis limits

        maxMat(boxNum,:) = [repmat(max(currData),[1,3])];
        minMat(boxNum,:) = [repmat(min(currData),[1,3])];
    
    elseif nnz(~isnan(currData)) <= 4

        boxMedian = median(currData,1,'omitnan'); % get median
        
        hold on;

        % plot median line
        line([xCoordinates(boxNum) - 0.25,xCoordinates(boxNum) + 0.25],...
             [boxMedian,boxMedian],...
             'Color',[medianColors{boxNum},medianAlpha],...
             'LineWidth',medianWidth,...
             'LineStyle',medianStyle,...
             'Tag','Median');

        % scatter outliers
        scatter(xCoordinates(boxNum) - jitterMat(:,boxNum),...
                currData,...
                outlierSize,...
                'MarkerFaceColor',outlierColors{boxNum},...
                'MarkerEdgeColor',outlierColors{boxNum},...
                'Marker',outlierStyle,...
                'MarkerFaceAlpha',outlierAlpha,...
                'MarkerEdgeAlpha',outlierAlpha,...
                'Tag','Outlier');

        %%%% store min/max from each category for setting axis limits

        maxMat(boxNum,:) = [repmat(max(currData),[1,3])];
        minMat(boxNum,:) = [repmat(min(currData),[1,3])];
        
    end
    
end

%%%% connect groups with lines

if (isnumeric(plotLines) || plotLines == true) && nGroups == nBoxes

    currData = reshape(inputData,[],nBoxes);

    if isnumeric(plotLines), lineIdx = plotLines; else, lineIdx = 1:nBoxes; end

    for connection = 1:numel(lineIdx) - 1

        for sample = find(~all(isnan(currData),2))'

            plot(xCoordinates(lineIdx(connection:connection + 1)) - jitterMat(sample,lineIdx(connection:connection + 1)),...
                 currData(sample,lineIdx(connection:connection + 1)),...
                 'Color',[lineColors{find(~all(isnan(currData),2)) == sample,connection},lineAlpha],...
                 'LineWidth',lineWidth,...
                 'LineStyle',lineStyle,...
                 'Tag','Line');
        
        end

    end

elseif (isnumeric(plotLines) || plotLines == true) && nGroups ~= nBoxes

    currData = reshape(inputData,[],nBoxes);

    if isnumeric(plotLines), lineIdx = plotLines; else, lineIdx = 1:nBoxes/nGroups; end

    connectionNum = 1;

    for group = 1:nGroups

        boxIdx = reshape(1:nBoxes,[],nGroups)'; boxIdx = boxIdx(:,lineIdx);

        for connection = 1:numel(boxIdx(group,:)) - 1

            for sample = find(~all(isnan(currData),2))'

                plot(xCoordinates(boxIdx(group,connection:connection + 1)) - jitterMat(sample,boxIdx(group,connection:connection + 1)),...
                     currData(sample,boxIdx(group,connection:connection + 1)),...
                     'Color',[lineColors{find(~all(isnan(currData),2)) == sample,connectionNum},lineAlpha],...
                     'LineWidth',lineWidth,...
                     'LineStyle',lineStyle,...
                     'Tag','Line');
            
            end

            connectionNum = connectionNum + 1;

        end

    end

end

%%%% plot individual data points

if (isnumeric(plotPoints) || plotPoints == true) && nGroups == nBoxes

    currData = reshape(inputData,[],nBoxes);

    if isnumeric(plotPoints), pointIdx = plotPoints; else, pointIdx = 1:nBoxes; end

    for box = pointIdx

        for sample = find(~all(isnan(currData),2))'

            scatter(xCoordinates(box) - jitterMat(sample,box),...
                    currData(sample,box),...
                    pointSize,...
                    'MarkerFaceColor',pointColors{find(~all(isnan(currData),2)) == sample,box},...
                    'MarkerEdgeColor',pointColors{find(~all(isnan(currData),2)) == sample,box},...
                    'Marker',pointStyle,...
                    'MarkerFaceAlpha',pointAlpha,...
                    'MarkerEdgeAlpha',pointAlpha,...
                    'Tag','Point');

        end

    end

elseif (isnumeric(plotPoints) || plotPoints == true) && nGroups ~= nBoxes

    currData = reshape(inputData,[],nBoxes);

    if isnumeric(plotPoints), pointIdx = plotPoints; else, pointIdx = 1:nBoxes/nGroups; end

    boxIdx = reshape(1:nBoxes,[],nGroups)'; boxIdx = boxIdx(:,pointIdx); 
    
    boxNum = 1;
    
    for group = 1:nGroups

        for box = 1:size(boxIdx,2)

            for sample = find(~all(isnan(currData),2))'

                scatter(xCoordinates(boxIdx(group,box)) - jitterMat(sample,boxIdx(group,box)),...
                        currData(sample,boxIdx(group,box)),...
                        pointSize,...
                        'MarkerFaceColor',pointColors{find(~all(isnan(currData),2)) == sample,boxNum},...
                        'MarkerEdgeColor',pointColors{find(~all(isnan(currData),2)) == sample,boxNum},...
                        'Marker',pointStyle,...
                        'MarkerFaceAlpha',pointAlpha,...
                        'MarkerEdgeAlpha',pointAlpha,...
                        'Tag','Point');

            end

            boxNum = boxNum + 1;
        
        end

    end

end

%%%% set x-axis limits

xlim([0.2*boxSpacing,xCoordinates(end) + (0.5*boxSpacing)]);

%%%% set x-axis tick labels

xLabPos = [];

if labelGroups == true

    for group = 1:(nBoxes/nGroups):nBoxes, xLabPos = horzcat(xLabPos,median(xCoordinates(group:group + (nBoxes/nGroups) - 1))); end

    set(gca,'xtick',xLabPos); set(gca,'xticklabel',boxLabels);

else
        
    set(gca,'xtick',xCoordinates); set(gca,'xticklabel',boxLabels);

end

%%%% set figure & axes appearance

set(gcf,'color','w');

set(gca,'box','off',...
        'XColor','k',...
        'YColor','k',...
        'TickDir','out',...
        'TickLength',[0.01,0.01],...
        'FontSize',13,...
        'LineWidth',1);

%%%% set y-axis limits

if ~all(isnan(maxMat),'all')

    yUpper = max(maxMat,[],'all'); yLower = min(minMat,[],'all');
    
    yExt = (yUpper - yLower)*0.2;
    
    ylim([yLower - yExt,yUpper + yExt]);

end

%%%% add legend if specified

if plotLegend == true && ~isempty(lgdLabels) && ~isempty(lgdColors)

    if lgdFontSize > 2, lgdLineHeight = (lgdFontSize - 2)*lgdLineHeight; end

    for lgdEntry = 1:numel(lgdColors)

        currData = inputData(inputLabels == lgdEntry,1); med = median(currData,1,'omitnan');

        % plot line with box color for legend
        line([xCoordinates(lgdEntry)-0.25 xCoordinates(lgdEntry)+0.25],...
             [med med],...
             'Color',lgdColors{lgdEntry},...
             'LineWidth',lgdLineHeight,...
             'Tag','Legend Line');
    
        % plot while line over line with box color for legend
        line([xCoordinates(lgdEntry)-0.25 xCoordinates(lgdEntry)+0.25],...
             [med med],...
             'Color',[1,1,1],...
             'LineWidth',lgdLineHeight,...
             'Tag','Legend Line Cover');

        set(gca,'Children',circshift(gca().Children,-2,1));

    end

    warning('off','MATLAB:handle_graphics:exceptions:SceneNode');

    lgdObject = legend(findobj(gca,'Tag','Legend Line'),...
                       strcat('\fontsize{',num2str(lgdFontSize),'}',lgdLabels),...
                       AutoUpdate = 'off',...
                       NumColumns = lgdColumns,...
                       Orientation = lgdOrientation,...
                       Box = lgdBox,...
                       Location = lgdLocation);

    lgdObject.ItemTokenSize = [30*lgdLineWidth,9];

    if any(lgdPosition), lgdObject.Position = lgdPosition; end

else
    
    lgdObject = [];

end

end
