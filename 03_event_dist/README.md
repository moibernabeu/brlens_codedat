# Event distance calculation

Run this for each splitted tree. There is an optimum value to parallelise using greasy, it is 4 cores. This configuration maximises the number of tasks that are executed per node.

```
./../src/event_dist.py -f ../01_get_trees/splitted/0076_0.txt -g data/0076_norm_groups.csv -o outputs -c 4
```

To join all the outputs:
```

```
