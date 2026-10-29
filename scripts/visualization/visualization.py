from visualization.boxplot_visualization import BoxplotVisualization
from visualization.comparison_visualization import ComparisonVisualization
from visualization.heatmap_visualization import HeatmapVisualization
from visualization.map_visualization import MapVisualization
from visualization.pixel_visualization import PixelVisualization
from visualization.timing_visualization import TimingVisualization

    
class PhenologyVisualization(
    HeatmapVisualization,
    BoxplotVisualization,
    PixelVisualization,
    MapVisualization,
    ComparisonVisualization,
    TimingVisualization,
):
    """Backward-compatible composite of all phenology graph families.

    New code can use a focused subclass such as :class:`MapVisualization`.
    The composite remains the default because the analysis pipeline creates
    one object and passes it to several graph-producing helpers.
    """

    graph_type = "all"
    graph_methods = frozenset().union(
        HeatmapVisualization.graph_methods,
        BoxplotVisualization.graph_methods,
        PixelVisualization.graph_methods,
        MapVisualization.graph_methods,
        ComparisonVisualization.graph_methods,
        TimingVisualization.graph_methods,
    )
            
