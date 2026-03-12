# Italian Quarterly National Accounts

A subset of the data used by Girolimetto et al. (2023) from the Italian
Quarterly National Accounts (output, income and expenditure sides)
spanning the period 2000:Q1-2019:Q4.

## Usage

``` r
# 21 time series of the Italian Quarterly National Accounts
itagdp

# 'agg_mat' and 'cons_mat' for the output side
outside

# 'agg_mat' and 'cons_mat' for the expenditure side
expside

# 'agg_mat' and 'cons_mat' for the income side
incside

# zero constraints matrix encompassing output, expenditure and income sides
gdpconsmat
```

## Format

`itagdp` is a \\(80 \times 21)\\ `ts` object, corresponding to 21 time
series of the Italian Quarterly National Accounts (2000:Q1-2019:Q4).

`outside`, `income` and `expenditure` are lists with two elements:

- `agg_mat` contains the \\(1 \times 2)\\, \\(2 \times 4)\\, or \\(6
  \times 8)\\ aggregation matrix according to output, income or
  expenditure side, respectively.

- `cons_mat` contains the \\(1 \times 3)\\, \\(2 \times 6)\\, or \\(6
  \times 14)\\ zero constraints matrix according to output, income or
  expenditure side, respectively.

`gdpconsmat` is the complete \\(9 \times 21)\\ zero constraints matrix
encompassing output, expenditure and income sides.

## Source

<https://ec.europa.eu/eurostat/web/national-accounts/>

## References

Girolimetto, D. and Di Fonzo, T. (2023), Point and probabilistic
forecast reconciliation for general linearly constrained multiple time
series, *Statistical Methods & Applications*, 33, 581-607.
[doi:10.1007/s10260-023-00738-6](https://doi.org/10.1007/s10260-023-00738-6)
.
