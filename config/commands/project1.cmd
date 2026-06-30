init 1000 1000

g_auto 10 990 10 990 2 5 2 5 -30 30 1000 1001 Random

generate

bin 132 All

wave 10

triangulate
voronoi
build_nav_graph 1 90.0 90.0
connect_to_graph 3 997 997 3 All

find_path_astar
save_metrics var/metrics/graph.csv

find_path_dekstra
save_metrics var/metrics/graph.csv

find_path_greedy
save_metrics var/metrics/graph.csv

end
