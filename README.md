# MarkovPatternTool
**Install via:**
```
pip3 install git+https://github.com/flomock/MarkovPatternTool.git@master
```



## How to use:

  **Example:**
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

**Options:**

command | what it does
  ------------- | -------------
-p, --path              |Path to fasta, csv, or directory with multiple files.  [required]
-n, --n_length   |Length of the resulting word.
-k, --k_length |Length of the analyzed tuples.
-f, --fasta             |Necessary if path contains fasta-file which should be used.
-r, --recursiv          |Uses also all subdirectories.
--log                   |Returns fold results with log scale. Easier interpretation of over-,under- occurrence.
--filter                |Filter out microsatellites;
--help                  |Show this message and exit.

<br><br>
<a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc-sa/4.0/">Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License</a>.
