# Seed to species distance calculation

Run this for each splitted tree. There is an optimum value to parallelise using greasy, it is 4 cores. This configuration maximises the number of tasks that are executed per node.

```
./../src/seq2seq_cladenorm.py -f ../01_get_trees/splitted/0005_19.txt -o outputs/ -p data/0005_norm_groups.csv -c 4
```

The outputs can be joined by using:
```
./../src/join_normalise.py
```
