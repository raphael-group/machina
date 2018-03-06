## Ovarian cancer

The following scripts generate the results for the McPherson *et al.* dataset.

* [`mcpherson_2016/run_pmh_sankoff.sh`](mcpherson_2016/run_pmh_sankoff.sh) enumerates all minimum migration histories
* [`mcpherson_2016/run_pmh.sh`](mcpherson_2016/run_pmh.sh) solves the PMH problem under various topological constraints
* [`mcpherson_2016/run_pmh_tr.sh`](mcpherson_2016/run_pmh_tr.sh) solves the PMH-TR problem under various topological constraints

## Breast cancer

The following scripts generate the results for the Hoadley *et al.* dataset.

* [`hoadley_2016/run_A1.sh`](hoadley_2016/run_A1.sh) solves the PMH-TI problem for patient A1. This is done by enumerating all mutation trees.
* [`hoadley_2016/run_A7.sh`](hoadley_2016/run_A7.sh) solves the PMH-TI problem for patient A7. This is done by enumerating all mutation trees and migration trees.

## Melonama

The following scripts generate the results for the Sanborn *et al.* dataset.

* [`sanborn_2015/run.sh`](sanborn_2015/run.sh) solves the PMH-TR problem for each melanoma patient.

## Prostate

The following scripts generate the results for the Gundem *et al.* dataset.

* [`gundem_2015/run.sh`](sanborn_2015/run.sh) solves the PMH-TR problem for each prostate patient.
* [`gundem_2015/run_A22_noCNA.sh`](gundem_2015/run_A22_noCNA.sh) solve the PMH-TI problem for patient A22 using only the reported copy-neutral SNVs.

## Simulations

The following scripts generate MACHINA results for the simulated instances.

* [`sims/machina/run.sh`](sims/machina/run.sh) solves the PMH-TI problem for each simulated instance.
* [`sims/machina/run_sankoff.sh`](sims/machina/run_sankoff.sh) enumerates all minimum migration histories of each simulated clone tree.

