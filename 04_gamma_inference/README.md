# Gamma inference

Here we stored the input and outputs of the Gamma inference process. In the `data` folder there are stored the filtered outputs of the distances calculations. In the `outputs` folder the MCMC outputs are stored in a 3 columns dataframe, the parameter (`param`), the `value` and the `event` or the destination species (`spto`).

To obtain these outputs `mcmc_sp_script.R` was executed for each species. The species specific scripts were generated automatically by changing the number of the species vector.

```
for i in {01..20}
do
  sed "s/sps\[1\]/sps\[$i\]/g" mcmc_sp_script.R > mcmc_yeast_${i}.R
done

for i in {01..38}
do
  sed "s/sps\[1\]/sps\[$i\]/g" mcmc_sp_script.R | sed 's/yedat/hudat/g' > mcmc_human_${i}.R
done
```
