"""
MPT command line.
"""

import click
import FractalMatrix

@click.command()
@click.option('--path','-p', required=True, help='Path to fasta, csv, or directory with multiple files.')
@click.option('--n_length','-n', default = 5, help='Length of the resulting word.')
@click.option('--k_length','-k', default = 2, help='Length of the analyzed tuples.')
@click.option('--fasta','-f', is_flag=True, help='Necessary if path contains fasta-file which should be used.')
@click.option('--recursiv','-r', is_flag=True, help='Uses also all subdirectories.')

# outputpath, save csv, plot csv-classic, plot csv-markov,
def cli(path,n_length,k_length,fasta,recursiv):
    '''
    Example:

    \b
    $ MPT -p /home/user/csv/arabidopsis_thaliana.csv
    \b
    use fasta-file
    $ MPT -f -p /home/user/fasta/arabidopsis_thaliana.fa
    \b
    use multiple fasta-files inclusive subdirectories
    $ MPT -f -r -p /home/user/fasta/

    '''
    if fasta:
        FractalMatrix.path_to_fastaFiles(path,recursiv)
    FractalMatrix.path_to_markovPatternAnalyse(path,n_length,k_length,recursiv)

