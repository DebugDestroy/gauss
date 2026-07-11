init 1000 700

g_auto 100 101 15 985 15 685 10 12 10 12 -20 20

generate

save_g

gnuplot results/visualizations/gnuplot.png
bmp_write results/visualizations/Pole.bmp Full

bin 5
bmp_write results/visualizations/Slice.bmp Binary

wave 10
PlotMetedata results/visualizations/Metadata.png

k_means 68
PlotKmeans results/visualizations/kmeans.png
k_means_kern 68 1
PlotKmeans results/visualizations/kmeans_with_kernel.png

triangulate
PlotDelaunay results/visualizations/Triangulation_Delone.png
voronoi
PlotVoronoi results/visualizations/Diagramma_Voronova.png
build_nav_graph 1 90 90
PlotGraph results/visualizations/Graph.png
connect_to_graph 2 2 998 698 All
PlotGraph results/visualizations/ConnectedGraph.png

astar_graph
save_metrics
PlotPath results/visualizations/AstarPathGraph.png
Plot3DPath results/visualizations/AstarPlot3DPathGraph.png

dekstra_graph
save_metrics
PlotPath results/visualizations/DekstraPathGraph.png
Plot3DPath results/visualizations/DekstraPlot3DPathGraph.png
greedy_graph
save_metrics
PlotPath results/visualizations/GreedyPathGraph.png
Plot3DPath results/visualizations/GreedyPlot3DPathGraph.png

grid 5
PlotGrid results/visualizations/Grid.png
build_nav_grid 0
PlotNavGrid results/visualizations/NavigationGrid.png
connect_to_grid 2 2 998 698
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

rrt 1000000 2 2 998 698 1 5 50 50 1 1 1 1 10 0.2
save_metrics
end
