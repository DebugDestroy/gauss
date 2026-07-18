init 50 50

g 1 1 2 2 1
g 1 15 2 2 1
g 15 1 2 2 1
g 20 11 2 2 1
g 10 17 2 2 1

generate
gnuplot results/visualizations/gnuplot.png
bmp_write results/visualizations/Pole.bmp Full

bin 1
bmp_write results/visualizations/Slice.bmp Binary

wave 0
PlotMetedata results/visualizations/Metadata.png

triangulate
PlotDelaunay results/visualizations/Triangulation_Delone.png
voronoi
PlotVoronoi results/visualizations/Diagramma_Voronova.png
build_nav_graph 1 90 90
PlotGraph results/visualizations/Graph.png
connect_to_graph 48 48 48 1 All
PlotGraph results/visualizations/ConnectedGraph.png

astar_graph
save_metrics
PlotPathDiscrete results/visualizations/AstarPathGraph.png
Plot3DPath results/visualizations/AstarPlot3DPathGraph.png

dekstra_graph
save_metrics
PlotPathDiscrete results/visualizations/DekstraPathGraph.png
Plot3DPath results/visualizations/DekstraPlot3DPathGraph.png
greedy_graph
save_metrics
PlotPathDiscrete results/visualizations/GreedyPathGraph.png
Plot3DPath results/visualizations/GreedyPlot3DPathGraph.png

grid 5
PlotGrid results/visualizations/Grid.png
build_nav_grid 0
PlotNavGrid results/visualizations/NavigationGrid.png
connect_to_grid 48 48 48 1
PlotNavGrid results/visualizations/ConnectedGrid.png

astar_grid
save_metrics
PlotGridPath results/visualizations/AstarPathGrid.png
dekstra_grid
save_metrics
PlotGridPath results/visualizations/DekstraPathGrid.png
greedy_grid
save_metrics
PlotGridPath results/visualizations/GreedyPathGrid.png

rrt 20 1000000 48 48 48 1 1 1 50 50 1 1 1 4 10 0.2
PlotPathContinuous results/visualizations/RRT.png
save_metrics
end
