# flowjo-python
Imports exported FlowJo CSV files and creates FlowJo gating tree from CSV headers. Calculates custom population frequencies, creates plots, and runs statistical tests

# Description 
Flow cytometry data is often high-dimensional, leading to large and complex gating schemes. After manual population gating scheme creation in FlowJo, the gated populations are easily exported as a CSV table of "Frequency of Parent" values. This Python package imports CSV tables and generates a hierarchal tree object (Bigtree) which describes the gating scheme created in FlowJo (dubbed "flowtree"). With this information, and "Event Count" totals, this package allows for population Counts and custom Frequencies to be calculated. 

Metadata included in the FlowJo CSVs (ie: Treatment Group) is saved to the flowtree. Additional data, such as MRTI or MSOT read-outs, can be appended to flowtree and merged on sample ID numbers. 

# In Development
- ggplot-style plots to compare immune populations across Treatment Group, or other metadata variables
- statistical tests (t-tests) across Treatment Group, or other metadata variables
- ggplot-style plots and correlation coefficient analysis with continuous variables, such as thermal exposure (MRTI) vs. T cell infiltration

# Dependencies 
[pandas 2.2.1](https://pandas.pydata.org/)
[bigtree 0.17.0](https://bigtree.readthedocs.io/en/latest/)
[plotnine 0.13.4](https://plotnine.org/)
[plotnine-prism](https://pwwang.github.io/plotnine-prism/)
[patchworklib](https://pypi.org/project/patchworklib/0.3.0/)
[seaborn](https://seaborn.pydata.org/)
