# MarkovPatternTool
Install via:
pip3 install git+https://github.com/flomock/MarkovPatternTool.git@master


How to use:

  Example:
```
  $ MPT -p /home/user/csv/arabidopsis_thaliana.csv
```
  use fasta-file
```
  $ MPT -f -p /home/user/fasta/arabidopsis_thaliana.fa
```
  use multiple fasta-files inclusive subdirectories
```
  $ MPT -f -r -p /home/user/fasta/ -k 3 -n 6
```
  run with multiple n or k values, eg k=0,1,2,3
```
  $ MPT -p /home/user/csv/arabidopsis_thaliana.csv -k 0-3
```  

Options:
  -p, --path TEXT         Path to fasta, csv, or directory with multiple
                          files.  [required]
  -n, --n_length INTEGER  Length of the resulting word.
  -k, --k_length INTEGER  Length of the analyzed tuples.
  -f, --fasta             Necessary if path contains fasta-file which should
                          be used.
  -r, --recursiv          Uses also all subdirectories.
  --help                  Show this message and exit.

