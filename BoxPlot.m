
classdef BoxPlot < matlab.graphics.chartcontainer.ChartContainer & ...
                   matlab.graphics.chartcontainer.mixin.Legend

    properties
        InputData
        InputLabels

        GroupSize
        GroupLabels % specify either labels themselves or boolean true
        
        BoxLabels
        BoxColors
        BoxWidth
        BoxAlpha
        BoxEdgeColors
        BoxEdgeWidth
        BoxEdgeStyle
        BoxEdgeAlpha
        BoxSpacing
        BoxCurvature
        
        MedianColors
        MedianWidth
        MedianStyle
        MedianAlpha
        
        WhiskerColors
        WhiskerWidth
        WhiskerStyle
        WhiskerAlpha
        WhiskerCapColors
        WhiskerCapWidth
        WhiskerCapStyle
        WhiskerCapAlpha
        WhiskerCapLength
        
        OutlierColors
        OutlierSize
        OutlierStyle
        OutlierAlpha
        OutlierJitter
        
        JitterWidth
        
        PointDisplay
        PointColors
        PointSize
        PointStyle
        PointAlpha
        PointJitter
        
        LineDisplay
        LineColors
        LineWidth
        LineStyle
        LineAlpha
        
        LegendLabels
        LegendColors
        LegendColumns
        LegendFontSize
        LegendLineHeight
        LegendLineWidth
        LegendOrientation
        LegendBox
        LegendLocation
        LegendPosition
    end
    properties (Access = private, Transient, NonCopyable)
        % OriginalBox
        % OriginalWhisker
        %
    end

    methods
        function obj = BoxPlot(input_data, varargin)
            % Check for at least one input.
            if nargin < 1
                error("Not enough inputs. Pass either:" + ...
                      "\n    1) A matrix of box data or" + ...
                      "\n    2) a vector with labels specified via " + ...
                               "the BoxLabels name-value argument.");
            end
            % Convert input_data into a name-value pair.
            args = {'InputData', input_data};
            % Combine args with user-provided name-value pairs.
            args = [args, varargin];
            % Call superclass constructor method.
            obj@matlab.graphics.chartcontainer.ChartContainer(args{:});
        end
    end
    methods (Access = protected)
        function setup(obj)
            ax = getAxes(obj);
        end
        function update(obj)
            ax = getAxes(obj);
        end
    end

end
