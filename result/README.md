## Ovarian cancer

The following scripts generate the results for the McPherson *et al.* dataset.

* [`mcpherson_2016/run_pmh_sankoff.sh`](mcpherson_2016/run_pmh_sankoff.sh) enumerates all minimum migration histories
* [`mcpherson_2016/run_pmh.sh`](mcpherson_2016/run_pmh.sh) solves the PMH problem under various topological constraints
* [`mcpherson_2016/run_pmh_pr.sh`](mcpherson_2016/run_pmh_pr.sh) solves the PMH-PR problem under various topological constraints

## Breast cancer

The following scripts generate the results for the Hoadley *et al.* dataset.

* [`hoadley_2016/run_A1.sh`](hoadley_2016/run_A1.sh) solves the PMH-CTI problem for patient A1. This is done by enumerating all migration trees. The minimum migration solution is obtained with migration tree 41.
* [`hoadley_2016/run_A7.sh`](hoadley_2016/run_A7.sh) solves the PMH-CTI problem for patient A7. This is done by enumerating all migration trees. The minimum migration solutions are obtained with migration trees 1040, 1045 and 1051.
