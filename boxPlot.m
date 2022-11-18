
function [xCoordinates, legendObject] = boxPlot(inputData, NameValueArgs)
%
% BOXPLOT Plot boxplot from vector or matrix input.
%
% Input Arguments:
%
%     inputData
%
% Name-Value Arguments (default):
%
%     inputLabels
%
%     groupSize
%     labelGroup
%
%     boxLabels
%     boxColors
%     boxAlpha
%     boxEdgeColors
%     boxEdgeWidth
%     boxEdgeStyle
%     boxEdgeAlpha
%     boxSpacing
%     boxCurvature
%
%     medianColors
%     medianWidth
%     medianStyle
%     medianAlpha
%
%     whiskerColors
%     whiskerWidth
%     whiskerStyle
%     whiskerAlpha
%
%     outlierColors
%     outlierSize
%     outlierStyle
%     outlierAlpha
%     outlierJitter
%
%     jitterWidth
%
%     pointDisplay
%     pointColors
%     pointSize
%     pointStyle
%     pointAlpha
%     pointJitter
%
%     lineDisplay
%     lineColors
%     lineWidth
%     lineStyle
%     lineAlpha
%
%     legendLabels
%     legendColors
%     legendColumns
%     legendFontSize
%     legendLineHeight
%     legendLineWidth
%     legendOrientation
%     legendBox
%     legendLocation
%     legendPosition
%
% Output Arguments:
%
%     xCoordinates
%
%     legendObject
%
% MIT License
% Copyright (c) 2022 Ryan Gorzek
% https://github.com/gorzek-ryan/matlab_viz/blob/main/LICENSE
% https://opensource.org/licenses/MIT
%

arguments
    
    inputData (:,:) {mustBeNumeric}
    
    NameValueArgs.inputLabels {mustBeVector, mustBeNumeric} = reshape(repmat(1:size(inputData,2),[size(inputData,1),1]),[],1)

    NameValueArgs.groupSize (1,1) {mustBeNumeric} = 1
    NameValueArgs.labelGroups (1,1) logical = false

    NameValueArgs.boxLabels {mustBeA(NameValueArgs.boxLabels,"cell")} = {}
    NameValueArgs.boxColors {mustBeA(NameValueArgs.boxColors,"cell")} = {}
    NameValueArgs.boxAlpha (1,1) {mustBeInRange(NameValueArgs.boxAlpha,0,1)} = 1
    NameValueArgs.boxEdgeColors {mustBeA(NameValueArgs.boxEdgeColors,"cell")} = {}
    NameValueArgs.boxEdgeWidth (1,1) {mustBeNumeric} = 1
    NameValueArgs.boxEdgeStyle (1,1) string = "-"
    NameValueArgs.boxEdgeAlpha (1,1) {mustBeInRange(NameValueArgs.boxEdgeAlpha,0,1)} = 1
    NameValueArgs.boxSpacing (1,1) {mustBeNumeric} = 1
    NameValueArgs.boxCurvature (1,2) double = [0,0]

    NameValueArgs.medianColors (1,:) {mustBeA(NameValueArgs.medianColors,"cell")} = {}
    NameValueArgs.medianWidth (1,1) {mustBeNumeric} = 2
    NameValueArgs.medianStyle (1,1) string = "-"
    NameValueArgs.medianAlpha (1,1) {mustBeInRange(NameValueArgs.medianAlpha,0,1)} = 1

    NameValueArgs.whiskerColors (1,:) {mustBeA(NameValueArgs.whiskerColors,"cell")} = {}
    NameValueArgs.whiskerWidth (1,1) {mustBeNumeric} = 1
    NameValueArgs.whiskerStyle (1,1) string = "-"
    NameValueArgs.whiskerAlpha (1,1) {mustBeInRange(NameValueArgs.whiskerAlpha,0,1)} = 1

    NameValueArgs.outlierColors (1,:) {mustBeA(NameValueArgs.outlierColors,"cell")} = {}
    NameValueArgs.outlierSize (1,1) {mustBeNumeric} = 30
    NameValueArgs.outlierStyle (1,1) string = "o"
    NameValueArgs.outlierAlpha (1,1) {mustBeInRange(NameValueArgs.outlierAlpha,0,1)} = 1
    NameValueArgs.outlierJitter {mustBeMember(NameValueArgs.outlierJitter,["none","rand","randn"])} = "none"

    NameValueArgs.jitterWidth (1,1) {mustBeNumeric} = 0.75

    NameValueArgs.pointDisplay {mustBeNumericOrLogical} = false
    NameValueArgs.pointColors {mustBeA(NameValueArgs.pointColors,"cell")} = {}
    NameValueArgs.pointSize (1,1) {mustBeNumeric} = 20
    NameValueArgs.pointStyle (1,1) string = "."
    NameValueArgs.pointAlpha (1,1) {mustBeInRange(NameValueArgs.pointAlpha,0,1)} = 1
    NameValueArgs.pointJitter {mustBeMember(NameValueArgs.pointJitter,["none","rand","randn"])} = "none"

    NameValueArgs.lineDisplay {mustBeNumericOrLogical} = false
    NameValueArgs.lineColors {mustBeA(NameValueArgs.lineColors,"cell")} = {}
    NameValueArgs.lineWidth (1,1) {mustBeNumeric} = 0.5
    NameValueArgs.lineStyle (1,1) string = "-"
    NameValueArgs.lineAlpha (1,1) {mustBeInRange(NameValueArgs.lineAlpha,0,1)} = 1

    NameValueArgs.legendLabels (1,:) {mustBeA(NameValueArgs.legendLabels,"cell")} = {}
    NameValueArgs.legendColors (1,:) {mustBeA(NameValueArgs.legendColors,"cell")} = {}
    NameValueArgs.legendColumns (1,1) {mustBeNumeric} = 1
    NameValueArgs.legendFontSize (1,1) {mustBeNumeric} = 12
    NameValueArgs.legendLineHeight (1,1) {mustBeNumeric} = 1
    NameValueArgs.legendLineWidth (1,1) {mustBeNumeric} = 1
    NameValueArgs.legendOrientation {mustBeMember(NameValueArgs.legendOrientation,["vertical","horizontal"])} = "horizontal"
    NameValueArgs.legendBox {mustBeMember(NameValueArgs.legendBox,["on","off"])} = "off"
    NameValueArgs.legendLocation (1,1) string = "northeast"
    NameValueArgs.legendPosition (1,4) = [0, 0, 0, 0]

end

% Reassign NameValueArgs structure to NVAs for short.
NVAs = NameValueArgs;

% Get number of boxes, groups, samples (including NaN), and missing values.
nBoxes = numel(unique(NVAs.inputLabels));
nGroups = nBoxes/NVAs.groupSize;
nSamples = size(inputData,1);
nMissing = nnz(all(isnan(inputData),2));

% Throw error if number of groups is not an integer.
if rem(nGroups,1) ~= 0
    error("Number of input boxes is not divisible by number of groups."); 
end

% Reshape input data into a column vector and get unique labels.
inputData = reshape(inputData,[],1);
uniqueLabels = sort(unique(NVAs.inputLabels),"ascend")';

% Set default box labels (numbered) if none are specified.
if isempty(NVAs.boxLabels) && ...
   NVAs.groupSize == 1

    NVAs.boxLabels = cellstr(string(1:nBoxes));

elseif isempty(NVAs.boxLabels) && ...
       NVAs.groupSize > 1 && ...
       NVAs.labelGroups == false

    NVAs.boxLabels = cellstr(string(1:NVAs.groupSize));

elseif isempty(NVAs.boxLabels) && ...
       NVAs.labelGroups == true
    
    NVAs.boxLabels = cellstr(string(1:nGroups));

end

% Set default box colors (gray) if none are specified, or replicate
% single or group-level color specification.
if isempty(NVAs.boxColors)
    NVAs.boxColors = repmat({[0.7,0.7,0.7]},[1,nBoxes]);
elseif numel(NVAs.boxColors) == 1
    NVAs.boxColors = repmat(NVAs.boxColors,[1,nBoxes]);
elseif numel(NVAs.boxColors) == NVAs.groupSize
    NVAs.boxColors = repmat(NVAs.boxColors,[1,nGroups]);
end

% Set default box edge colors (black) if none are specified, or replicate
% single or group-level color specification.
if isempty(NVAs.boxEdgeColors)
    NVAs.boxEdgeColors = repmat({[0.0,0.0,0.0]},[1,nBoxes]);
elseif numel(NVAs.boxEdgeColors) == 1
    NVAs.boxEdgeColors = repmat(NVAs.boxEdgeColors,[1,nBoxes]);
elseif numel(NVAs.boxEdgeColors) == NVAs.groupSize
    NVAs.boxEdgeColors = repmat(NVAs.boxEdgeColors,[1,nGroups]);
end

% Set default whisker colors (black) if none are specified, or replicate
% single or group-level color specification.
if isempty(NVAs.whiskerColors)
    NVAs.whiskerColors = repmat({[0.0,0.0,0.0]},[1,nBoxes]);
elseif numel(NVAs.whiskerColors) == 1
    NVAs.whiskerColors = repmat(NVAs.whiskerColors,[1,nBoxes]);
elseif numel(NVAs.whiskerColors) == NVAs.groupSize
    NVAs.whiskerColors = repmat(NVAs.whiskerColors,[1,nGroups]);
end

% Set default median colors (black) if none are specified, or replicate
% single or group-level color specification.
if isempty(NVAs.medianColors)
    NVAs.medianColors = repmat({[0.0,0.0,0.0]},[1,nBoxes]);
elseif numel(NVAs.medianColors) == 1
    NVAs.medianColors = repmat(NVAs.medianColors,[1,nBoxes]);
elseif numel(NVAs.medianColors) == NVAs.groupSize
    NVAs.medianColors = repmat(NVAs.medianColors,[1,nGroups]);
end

% Set default outlier colors (matched to box colors) if none are specified,
% or replicate single or group-level color specification.
if isempty(NVAs.outlierColors)
    NVAs.outlierColors = NVAs.boxColors;
elseif numel(NVAs.outlierColors) == 1
    NVAs.outlierColors = repmat(NVAs.outlierColors,[1,nBoxes]);
elseif numel(NVAs.outlierColors) == NVAs.groupSize
    NVAs.outlierColors = repmat(NVAs.outlierColors,[1,nGroups]);
end

% Set default point colors (black) if not specified.
if isempty(NVAs.pointColors)
    NVAs.pointColors = repmat({[0.0,0.0,0.0]},[nSamples,nBoxes]);
elseif numel(NVAs.pointColors) == 1
    NVAs.pointColors = repmat(NVAs.pointColors,[nSamples,nBoxes]);
elseif isvector(NVAs.pointColors) && ...
       any(size(NVAs.pointColors) == NVAs.groupSize)
    NVAs.pointColors = repmat(NVAs.pointColors,[nSamples,nGroups]);
elseif isvector(NVAs.pointColors) && ...
       any(size(NVAs.pointColors) == nBoxes)
    NVAs.pointColors = repmat(NVAs.pointColors,[nSamples,1]);
elseif all(ismember(size(NVAs.pointColors),[nSamples - nMissing,NVAs.groupSize]))
    NVAs.pointColors = repmat(NVAs.pointColors,[1,nGroups]);
end

% Check default color stream specification for point colors.
point_defaultIdx = cell2mat(cellfun(@(x) strcmp(x,"default"),NVAs.pointColors,"UniformOutput",false));
if any(point_defaultIdx,"all")
    defaultColors = repmat(num2cell(colororder,2),[ceil(size(NVAs.pointColors,1)/7), ...
                                                   size(NVAs.pointColors,2)]);
    defaultColors = defaultColors(1:size(NVAs.pointColors,1),:);
    NVAs.pointColors(point_defaultIdx) = defaultColors(point_defaultIdx);
end

% Set default line colors (black) if not specified.
nLines = nGroups*(NVAs.groupSize - 1);
if isempty(NVAs.lineColors)
    NVAs.lineColors = repmat({[0.0,0.0,0.0]},[nSamples,nLines]);
elseif numel(NVAs.lineColors) == 1
    NVAs.lineColors = repmat(NVAs.lineColors,[nSamples,nLines]);
elseif isvector(NVAs.lineColors) && ...
       any(size(NVAs.lineColors) == nGroups)
    NVAs.lineColors = NVAs.lineColors(kron(1:nGroups,ones(nSamples,groupSize-1)));
elseif isvector(NVAs.lineColors) && ...
       any(size(NVAs.lineColors) == nLines)
    NVAs.lineColors = repmat(NVAs.lineColors,[nSamples,1]);
elseif all(ismember(size(NVAs.lineColors),[nSamples - nMissing,nGroups]))
    NVAs.lineColors = NVAs.lineColors(:,kron(1:nGroups,ones(1,groupSize-1)));
end

% Check default color stream specification for line colors.
line_defaultIdx = cell2mat(cellfun(@(x) strcmp(x,"default"),NVAs.lineColors,"UniformOutput",false));
if any(line_defaultIdx,"all")
    defaultColors = repmat(num2cell(colororder,2),[ceil(size(NVAs.lineColors,1)/7), ...
                                                   size(NVAs.lineColors,2)]);
    defaultColors = defaultColors(1:size(NVAs.lineColors,1),:);
    NVAs.lineColors(line_defaultIdx) = defaultColors(line_defaultIdx);
end

% Check boxLabels against nGroups to warn user about specifying labelGroups.
if NVAs.labelGroups == true && ...
   numel(NVAs.boxLabels) > nGroups
    error(["Number of box labels exceeds number of groups, " ...
           "did you mean to specify labelGroups = false?"]);
elseif NVAs.labelGroups == true && ...
       numel(NVAs.boxLabels) < nGroups
    error("Insufficient number of box labels for number of groups.");
elseif NVAs.labelGroups == false && ...
       NVAs.groupSize ~= 1 && ...
       numel(NVAs.boxLabels) == nGroups && ...
       numel(NVAs.boxColors) == NVAs.groupSize
    error(["Insufficient number of box labels for number of groups, " ...
           "did you mean to specify labelGroups = true?"]);
end

% Generate x-coordinates for plotting boxes.
if nGroups == 1
    xCoordinates = (0.70 : 0.65 : 0.70 + (0.65*nBoxes) - 0.65).*boxSpacing;
else
    xCoordinates = 0.70 : 0.55 : 0.70 + (0.55*(nBoxes/nGroups) - 0.55); 
    initCoors = 0.70 : 0.55 : 0.70 + (0.55*(nBoxes/nGroups) - 0.55);
    for grp = 2:nGroups
        xCoordinates = horzcat(xCoordinates, initCoors + xCoordinates(end) + 0.4);
    end
    xCoordinates = xCoordinates.*NVAs.boxSpacing;
end

% Generate matrix of jitter if specified.
if strcmp(NVAs.outlierJitter,"rand") || ...
   strcmp(NVAs.pointJitter,"rand")

    jitterMat = rand(size(reshape(inputData,[],nBoxes))).*(0.5*NVAs.jitterWidth);
    jitterMat = jitterMat - mean(jitterMat,1);

elseif strcmp(NVAs.outlierJitter,"randn") || ...
       strcmp(NVAs.pointJitter,"randn")

    randnMat = randn(size(reshape(inputData,[],nBoxes)));
    jitterMat = (randnMat./max(randnMat,[],1)).*(0.5*NVAs.jitterWidth);
    jitterMat = jitterMat - mean(jitterMat,1);

elseif strcmp(NVAs.outlierJitter,"none") && ...
       strcmp(NVAs.pointJitter,"none")

    jitterMat = zeros(size(reshape(inputData,[],nBoxes)));

end

% Intialize matrices for storing max/min values to set axis limits.
maxMat = zeros(nBoxes,3); minMat = zeros(nBoxes,3);

% Plot boxes.
for box = uniqueLabels

    % Get current box location and data.
    boxNum = find(uniqueLabels == box);
    currData = inputData(NVAs.inputLabels == box,1);
    
    % Plot box if there are at least 4 data points.
    if nnz(~isnan(currData)) > 4

        boxMedian = median(currData,1,"omitnan");

        upperQuantile = quantile(currData, 0.75);
        lowerQuantile = quantile(currData, 0.25);

        maxWhisker = upperQuantile + 1.5*(upperQuantile - lowerQuantile); 
        upperWhisker = max(currData(currData < maxWhisker & currData >= upperQuantile));
        if isempty(upperWhisker), upperWhisker = maxWhisker; end

        minWhisker = lowerQuantile - 1.5*(upperQuantile - lowerQuantile);
        lowerWhisker = min(currData(currData > minWhisker & currData <= lowerQuantile));
        if isempty(lowerWhisker), lowerWhisker = minWhisker; end
        
        hold on;

        % Plot box with bounds at quartiles.
        rectangle("Position",  [xCoordinates(boxNum)-0.25, lowerQuantile, 0.5, upperQuantile-lowerQuantile],...
                  "FaceColor", [NVAs.boxColors{boxNum}, NVAs.boxAlpha],...
                  "EdgeColor", [NVAs.boxEdgeColors{boxNum}, NVAs.boxEdgeAlpha],...
                  "LineWidth",  NVAs.boxEdgeWidth,...
                  "LineStyle",  NVAs.boxEdgeStyle,...
                  "Curvature",  NVAs.boxCurvature,...
                  "Tag",       "Box");

        % Plot median line.
        line([xCoordinates(boxNum)-0.25, xCoordinates(boxNum)+0.25],...
             [boxMedian, boxMedian],...
             "Color",    [NVAs.medianColors{boxNum}, NVAs.medianAlpha],...
             "LineWidth", NVAs.medianWidth,...
             "LineStyle", NVAs.medianStyle,...
             "Tag",      "Median");

        % Plot lower whisker.
        line([xCoordinates(boxNum), xCoordinates(boxNum)],...
             [lowerQuantile, lowerWhisker],...
             "Color",    [NVAs.whiskerColors{boxNum}, NVAs.whiskerAlpha],...
             "LineWidth", NVAs.whiskerWidth,...
             "LineStyle", NVAs.whiskerStyle,...
             "Tag",      "Lower Whisker");

        % Plot upper whisker.
        line([xCoordinates(boxNum), xCoordinates(boxNum)],...
             [upperQuantile, upperWhisker],...
             "Color",    [NVAs.whiskerColors{boxNum}, NVAs.whiskerAlpha],...
             "LineWidth", NVAs.whiskerWidth,...
             "LineStyle", NVAs.whiskerStyle,...
             "Tag",      "Upper Whisker");

        % Plot lower whisker bar.
        line([xCoordinates(boxNum)-0.1, xCoordinates(boxNum)+0.1],...
             [lowerWhisker, lowerWhisker],...
             "Color",    [NVAs.whiskerColors{boxNum}, NVAs.whiskerAlpha],...
             "LineWidth", NVAs.whiskerWidth,...
             "LineStyle", NVAs.whiskerStyle,...
             "Tag",      "Lower Whisker Bar");

        % Plot upper whisker bar.
        line([xCoordinates(boxNum)-0.1, xCoordinates(boxNum)+0.1],...
             [upperWhisker, upperWhisker],...
             "Color",    [NVAs.whiskerColors{boxNum}, NVAs.whiskerAlpha],...
             "LineWidth", NVAs.whiskerWidth,...
             "LineStyle", NVAs.whiskerStyle,...
             "Tag",      "Upper Whisker Bar");

        % Plot upper outliers.
        if any(currData > upperWhisker)
            scatter(xCoordinates(boxNum) - jitterMat(currData > upperWhisker,boxNum),...
                    currData(currData > upperWhisker),...
                    NVAs.outlierSize,...
                    "MarkerFaceColor", NVAs.outlierColors{boxNum},...
                    "MarkerEdgeColor", NVAs.outlierColors{boxNum},...
                    "Marker",          NVAs.outlierStyle,...
                    "MarkerFaceAlpha", NVAs.outlierAlpha,...
                    "MarkerEdgeAlpha", NVAs.outlierAlpha,...
                    "Tag",            "Outlier");
            upperOutliers = currData(currData > upperWhisker);
        else
            upperOutliers = nan;
        end
        
        % Plot lower outliers.
        if any(currData < lowerWhisker)

            scatter(xCoordinates(boxNum) - jitterMat(currData < lowerWhisker,boxNum),...
                    currData(currData < lowerWhisker),...
                    NVAs.outlierSize,...
                    "MarkerFaceColor", NVAs.outlierColors{boxNum},...
                    "MarkerEdgeColor", NVAs.outlierColors{boxNum},...
                    "Marker",          NVAs.outlierStyle,...
                    "MarkerFaceAlpha", NVAs.outlierAlpha,...
                    "MarkerEdgeAlpha", NVAs.outlierAlpha,...
                    "Tag",            "Outlier");

            lowerOutliers = currData(currData < lowerWhisker);
        else
            lowerOutliers = nan;
        end

        % Store min/max from each box to set axis limits.
        maxMat(boxNum,:) = [upperQuantile, upperWhisker, max(upperOutliers)];
        minMat(boxNum,:) = [lowerQuantile, lowerWhisker, min(lowerOutliers)];
        
    elseif nnz(~isnan(currData)) == 0
        
        % Store min/max from each box to set axis limits.
        maxMat(boxNum,:) = [repmat(max(currData),[1,3])];
        minMat(boxNum,:) = [repmat(min(currData),[1,3])];
    
    elseif nnz(~isnan(currData)) <= 4

        boxMedian = median(currData,1,"omitnan");
        
        hold on;

        % Plot median line.
        line([xCoordinates(boxNum)-0.25, xCoordinates(boxNum)+0.25],...
             [boxMedian, boxMedian],...
             "Color",    [NVAs.medianColors{boxNum}, NVAs.medianAlpha],...
             "LineWidth", NVAs.medianWidth,...
             "LineStyle", NVAs.medianStyle,...
             "Tag",      "Median");

        % Scatter outliers.
        scatter(xCoordinates(boxNum) - jitterMat(:,boxNum),...
                currData,...
                outlierSize,...
                "MarkerFaceColor", NVAs.outlierColors{boxNum},...
                "MarkerEdgeColor", NVAs.outlierColors{boxNum},...
                "Marker",          NVAs.outlierStyle,...
                "MarkerFaceAlpha", NVAs.outlierAlpha,...
                "MarkerEdgeAlpha", NVAs.outlierAlpha,...
                "Tag",            "Outlier");

        % Store min/max from each box to set axis limits.
        maxMat(boxNum,:) = [repmat(max(currData), [1,3])];
        minMat(boxNum,:) = [repmat(min(currData), [1,3])];
        
    end
    
end

% Connect groups with lines if specified.
if (isnumeric(NVAs.plotLines) || NVAs.plotLines == true) && ...
   nGroups == nBoxes

    currData = reshape(inputData,[],nBoxes);

    if isnumeric(NVAs.plotLines)
        lineIdx = NVAs.plotLines; 
    else
        lineIdx = 1:nBoxes; 
    end

    for connection = 1:numel(lineIdx)-1
        for sample = find(~all(isnan(currData),2))'

            plot(xCoordinates(lineIdx(connection:connection+1)) - jitterMat(sample,lineIdx(connection:connection+1)),...
                 currData(sample,lineIdx(connection:connection+1)),...
                 "Color",    [NVAs.lineColors{find(~all(isnan(currData),2)) == sample,connection}, NVAs.lineAlpha],...
                 "LineWidth", NVAs.lineWidth,...
                 "LineStyle", NVAs.lineStyle,...
                 "Tag",      "Line");

        end
    end

elseif (isnumeric(NVAs.plotLines) || NVAs.plotLines == true) && ...
       nGroups ~= nBoxes

    currData = reshape(inputData,[],nBoxes);

    if isnumeric(NVAs.plotLines)
        lineIdx = NVAs.plotLines;
    else
        lineIdx = 1:nBoxes/nGroups;
    end

    connectionNum = 1;
    for group = 1:nGroups

        boxIdx = reshape(1:nBoxes,[],nGroups)'; 
        boxIdx = boxIdx(:,lineIdx);

        for connection = 1:numel(boxIdx(group,:))-1
            for sample = find(~all(isnan(currData),2))'

                plot(xCoordinates(boxIdx(group,connection:connection + 1)) - jitterMat(sample,boxIdx(group,connection:connection + 1)),...
                     currData(sample,boxIdx(group,connection:connection + 1)),...
                     "Color",    [NVAs.lineColors{find(~all(isnan(currData),2)) == sample,connectionNum}, NVAs.lineAlpha],...
                     "LineWidth", NVAs.lineWidth,...
                     "LineStyle", NVAs.lineStyle,...
                     "Tag",      "Line");
            
            end
            connectionNum = connectionNum + 1;
        end

    end

end

% Plot individual data points if specified.
if (isnumeric(NVAs.plotPoints) || NVAs.plotPoints == true)

    currData = reshape(inputData,[],nBoxes);

    if isnumeric(NVAs.plotPoints)
        pointIdx = NVAs.plotPoints;
    else
        pointIdx = 1:nBoxes;
    end

    for sample = find(~all(isnan(currData),2))'

        scatter(xCoordinates(pointIdx) - jitterMat(sample,pointIdx),...
                currData(sample,pointIdx),...
                NVAs.pointSize,...
                vertcat(NVAs.pointColors{find(~all(isnan(currData),2)) == sample,pointIdx}),...
                "Marker",          NVAs.pointStyle,...
                "MarkerFaceAlpha", NVAs.pointAlpha,...
                "MarkerEdgeAlpha", NVAs.pointAlpha,...
                "Tag",            "Point");

    end

end

% Set x-axis limits.
xlim([0.2*NVAs.boxSpacing, xCoordinates(end)+(0.5*NVAs.boxSpacing)]);

% Set x-axis tick labels.
xLabPos = [];
if NVAs.labelGroups == true
    for group = 1:(nBoxes/nGroups):nBoxes
        xLabPos = horzcat(xLabPos, median(xCoordinates(group : group+(nBoxes/nGroups)-1)));
    end
    set(gca,"xtick",      xLabPos,...
            "xticklabel", NVAs.boxLabels);
else 
    set(gca,"xtick",      xCoordinates,...
            "xticklabel", NVAs.boxLabels);
end

% Set figure and axes appearance.
set(gcf,"color",      "w");
set(gca,"box",        "off",...
        "XColor",     "k",...
        "YColor",     "k",...
        "TickDir",    "out",...
        "TickLength", [0.01,0.01],...
        "FontSize",    13,...
        "LineWidth",   1);

% Set y-axis limits.
if ~all(isnan(maxMat),"all")
    yUpper = max(maxMat,[],"all"); yLower = min(minMat,[],"all");
    yExt = (yUpper - yLower)*0.2;
    ylim([yLower-yExt, yUpper+yExt]);
end

% Plot legend if specified.
if ~isempty(NVAs.legendLabels) && ...
   ~isempty(NVAs.legendColors)

    if NVAs.legendFontSize > 2
        NVAs.legendLineHeight = (NVAs.legendFontSize - 2) * NVAs.legendLineHeight;
    end

    for legendEntry = 1:numel(NVAs.legendColors)

        currData = inputData(NVAs.inputLabels == legendEntry,1);
        med = median(currData,1,"omitnan");

        % Plot line with box color for legend.
        line([xCoordinates(legendEntry)-0.25, xCoordinates(legendEntry)+0.25],...
             [med, med],...
             "Color",     NVAs.legendColors{legendEntry},...
             "LineWidth", NVAs.legendLineHeight,...
             "Tag",      "Legend Line");
    
        % Plot white line over line with box color for legend.
        line([xCoordinates(legendEntry)-0.25, xCoordinates(legendEntry)+0.25],...
             [med, med],...
             "Color",    [1.0,1.0,1.0],...
             "LineWidth", NVAs.legendLineHeight,...
             "Tag",      "Legend Line Cover");

        set(gca, "Children",circshift(gca().Children,-2,1));

    end

    warning("off","MATLAB:handle_graphics:exceptions:SceneNode");

    legendObject = legend(findobj(gca, "Tag","Legend Line"),...
                       strcat('\fontsize{',num2str(NVAs.legendFontSize),'}',NVAs.legendLabels),...
                       "AutoUpdate", "off",...
                       "NumColumns",  NVAs.legendColumns,...
                       "Orientation", NVAs.legendOrientation,...
                       "Box",         NVAs.legendBox,...
                       "Location",    NVAs.legendLocation);

    legendObject.ItemTokenSize = [30 * NVAs.legendLineWidth, 9];

    if any(NVAs.legendPosition)
        legendObject.Position = NVAs.legendPosition;
    end

else
    
    legendObject = [];

end

end
