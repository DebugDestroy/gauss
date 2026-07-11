init 1000 1000

g_auto 900 901 10 990 10 990 3 6 3 6 -20 20

generate

bin 5

wave 10

triangulate
voronoi
build_nav_graph 1 30.0 30.0
connect_to_graph 3 997 997 3 All

astar_graph
save_metrics var/metrics/graph.csv

dekstra_graph
save_metrics var/metrics/graph.csv

greedy_graph
save_metrics var/metrics/graph.csv

end
