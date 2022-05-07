# cluster anaysis

```
usage: Cluster.py [-h] [--output OUTPUT] [--prefix PREFIX] [--score SCORE] [--test] config

Clustering pdb files within a path

positional arguments:
  config                Config file for clustering in JSON format.

optional arguments:
  -h, --help            show this help message and exit
  --output OUTPUT, -o OUTPUT
                        The output directory.
  --prefix PREFIX, -p PREFIX
                        The prefix of collected pdb files. (model/multimer)
  --score SCORE, -s SCORE
                        The score type used to rank classes. (ptm/plddt)
  --test                Wether to use test mode.
```

The ouptdir contains all ranked class folders, `cluster_result.csv` for human reading and `cluster_result.json` which has all clustering information.

In `config.json`:

```
{
    "inference_dir": "/mnt/vepfs/casp15/casp15_inference/T1106s1",
    "gromacs":{
        "cluster_cutoff": "0.2",
        "cluster_method": "gromos",
        "cluster_group" : "protein-h",
        "clean": true
    },
    "comment.0": "Units for kapps is kJ/(mol*A^2)"
}
```

`inference_dir` specifies the folder containing predicted pdf files.



## Usage Example

* To cluster multimers of T1109 and rank them by ptm, use JSON as :

```
{
    "inference_dir": "/mnt/vepfs/casp15/casp15_inference/T1109",
    "gromacs":{
        "cluster_cutoff": "0.15",
        "cluster_method": "gromos",
        "cluster_group" : "protein-h",
        "forcefield": "amber99sb",
        "clean": true
    },
}
```

use command:
```
python Cluster.py config_T1109.json -o T1109_cluster_result -p multimer -s ptm
```