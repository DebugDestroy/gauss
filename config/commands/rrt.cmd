init 1000 1000

g_auto 900 901 10 990 10 990 3 6 3 6 -20 20

generate

rrt 1000000 3 997 997 3 1 5 30 30 2 2 2 6 10 0.2

save_metrics var/metrics/rrt.csv

end
