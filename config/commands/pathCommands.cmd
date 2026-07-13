init 1000 1000

g_auto 900 901 10 990 10 990 3 6 3 6 -20 20

generate

bin 5

wave 10

triangulate
voronoi
build_nav_graph 1 50 50
connect_to_graph 2 2 998 998 All

astar_graph
shortcut_discrete
spline_discrete 20 1 5 50 50 2 2 2
save_metrics var/metrics/pathMetrics.csv
PlotPathContinuous results/visualizations/AstarPathGraphSpline.png
dekstra_graph
shortcut_discrete
spline_discrete 20 1 5 50 50 2 2 2
save_metrics var/metrics/pathMetrics.csv
PlotPathContinuous results/visualizations/DekstraPathGraphSpline.png
greedy_graph
shortcut_discrete
spline_discrete 20 1 5 50 50 2 2 2
save_metrics var/metrics/pathMetrics.csv
PlotPathContinuous results/visualizations/GreedyPathGraphSpline.png


grid 5
build_nav_grid 5
connect_to_grid 2 2 998 998

astar_grid
shortcut_discrete
spline_discrete 20 1 5 50 50 2 2 2
save_metrics var/metrics/pathMetrics.csv
PlotPathContinuous results/visualizations/AstarPathGridSpline.png
dekstra_grid
shortcut_discrete
spline_discrete 20 1 5 50 50 2 2 2
save_metrics var/metrics/pathMetrics.csv
PlotPathContinuous results/visualizations/DekstraPathGridSpline.png
greedy_grid
shortcut_discrete
spline_discrete 20 1 5 50 50 2 2 2
save_metrics var/metrics/pathMetrics.csv
PlotPathContinuous results/visualizations/GreedyPathGridSpline.png


rrt 20000 2 2 998 998 1 5 50 50 2 2 2 40 10 0.2
PlotPathContinuous results/visualizations/RRT.png
shortcut_continuous
PlotPathContinuous results/visualizations/RRTShortcut.png
spline_continuous 20
PlotPathContinuous results/visualizations/RRTSpline.png
save_metrics var/metrics/pathMetrics.csv
rrt_star 20000 2 2 998 998 1 5 50 50 2 2 2 40 80 700 10 0.00005
PlotPathContinuous results/visualizations/RRT*.png
shortcut_continuous
PlotPathContinuous results/visualizations/RRTStarShortcut.png
spline_continuous 20
PlotPathContinuous results/visualizations/RRTStarSpline.png
save_metrics var/metrics/pathMetrics.csv
end
