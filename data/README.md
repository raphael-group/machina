We consider two datasets:

* Metastatic ovarian cancer [[1]](#ref1)
* Metastatic breast cancer [[2]](#ref2)

## Ovarian cancer

McPherson *et al.* [[1]](#ref1) analyzed the clonal history of seven metastatic ovarian cancer patients. The authors used [PyClone](https://bitbucket.org/aroth85/pyclone/wiki/Home) to infer mutation clusters, and used these clusters to infer clone trees. Subsequently, the authors used the Sankoff algorithm to infer a single vertex labeling that minimizes the number of migrations.

In [`data/mcpherson_2016/`](data/mcpherson_2016/) we provide the reported clone trees as well as the reported vertex labeled. For patient1 for instance we have the following files:

* [`patient1.tree`](data/mcpherson_2016/patient1.tree) is the reported clone tree
* [`patient1.labeling`](data/mcpherson_2016/patient1.labeling) is the reported leaf labeling
* [`patient1.reported.labeling`](data/mcpherson_2016/patient1.reported.labeling) is the reported vertex labeling.

## Breast cancer

Hoadley *et al.* [[2]](#ref2) analyzed two metastatic breast cancer patients ([A1](data/hoadley_2016/A1) and [A7](data/hoadley_2016/A7)). The authors focused on only the copy-neutral single-nucleotide variants, which they clustered using [SciClone](https://github.com/genome/sciclone).

Using [SPRUCE](https://github.com/raphael-group/spruce) and the provided read counts, we infer four mutation trees for A1 and two mutation trees for A7. See [`run_SPRUCE.sh`](data/hoadley_2016/run_SPRUCE) for the exact commands. The ipython notebooks for inferring the 99.99% confidence intervals on mutation cluster frequencies are present in the `raw` subdirectories.

## References
<a name="ref1"></a>
[1] McPherson *et al.* Divergent modes of clonal spread and intraperitoneal mixing in high-grade serous ovarian cancer. [*Nature Genetics*, May 2016.](http://www.nature.com/ng/journal/v48/n7/abs/ng.3573.html)
<a name="ref2"></a>
[2] Hoadley *et al.* Tumor Evolution in Two Patients with Basal-like Breast Cancer: A Retrospective Genomics Study of Multiple Metastases. [*PLOS Med*, 13(12):e1002174, December 2016.](http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1002174)

