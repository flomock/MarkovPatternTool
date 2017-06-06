# MarkovPatternTool
Install via:
pip install git+https://github.com/flomock/MarkovPatternTool.git@master


How to use:

  Example:

  $ MPT -p /home/user/csv/arabidopsis_thaliana.csv
  
  use fasta-file
  $ MPT -f -p /home/user/fasta/arabidopsis_thaliana.fa
  
  use multiple fasta-files inclusive subdirectories
  $ MPT -f -r -p /home/user/fasta/

Options:
  -p, --path TEXT         Path to fasta, csv, or directory with multiple
                          files.  [required]
  -n, --n_length INTEGER  Length of the resulting word.
  -k, --k_length INTEGER  Length of the analyzed tuples.
  -f, --fasta             Necessary if path contains fasta-file which should
                          be used.
  -r, --recursiv          Uses also all subdirectories.
  --help                  Show this message and exit.

