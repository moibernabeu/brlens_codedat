# Getting the trees

`get_trees.py` allows to get the trees from PhylomeDB in parallel.

```
./../src/get_trees.py -f data/phylome_list.txt -w outputs -t
```

Posteriorly, to get equal size files, the script `spilt_trees.py` should be run. It returns a set of multiple tree files.

```
./../src/split_trees.py
```
