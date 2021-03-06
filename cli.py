"""
MPT command line.
"""

import click
import FractalMatrix
import os


@click.command()
@click.option('--path', '-p', required=True, help='Path to fasta, csv, or directory with multiple files.')
@click.option('--threads', '-t', default=str(0), help='Number of threads to use, default equal to number of CPU cores.')
@click.option('--n_length', '-n', default=str(5), help='Length of the resulting word.')
@click.option('--k_length', '-k', default=str(2), help='Length of the analyzed tuples.')
@click.option('--fasta', '-f', is_flag=True, help='Necessary if path contains fasta-file which should be used.')
@click.option('--recursiv', '-r', is_flag=True, help='Uses also all subdirectories.')
@click.option('--log', is_flag=True,
              help='Returns fold results with log scale. Easier interpretation of over-,under- occurrence.')
@click.option('--filter', is_flag=True, help='Filter out microsatellites')
@click.version_option(version=0.25, prog_name="Markov Pattern Tool")
# @click.option('--sat_length', default=str(5), help='Maximum length microsatellite tuple.')
# @click.option('--sat_count', default=str(10), help='Minimum number repetitions microsatellite tuple.')

def cli(path, threads, n_length, k_length, fasta, recursiv, log, filter):
    '''
    Example:

    \b
    $ MPT -p /home/user/csv/arabidopsis_thaliana.csv
    \b
    use fasta-file
    $ MPT -f -p /home/user/fasta/arabidopsis_thaliana.fa
    \b
    use multiple fasta-files inclusive subdirectories
    $ MPT -f -r -p /home/user/fasta/ -k 3 -n 6
    \b
    run with multiple n or k values, eg k=0,1,2,3
    $ MPT -p /home/user/csv/arabidopsis_thaliana.csv -k 0-3
    \b
    output results with logarithmic scale
    $ MPT --log -f -r -p /home/user/fasta/ -k 3 -n 6
    '''
    if ('-' in n_length):
        n_start = int(str(n_length).split('-')[0])
        n_stop = int(str(n_length).split('-')[1])
    else:
        n_start = n_stop = int(n_length)

    if ('-' in k_length):
        k_start = int(str(k_length).split('-')[0])
        k_stop = int(str(k_length).split('-')[1])
    else:
        k_start = k_stop = int(k_length)

    if fasta:
        if (n_stop > 8):
            FractalMatrix.path_to_fastaFiles(path, recursiv, n=n_stop, filter_mikroSats=filter, threads=threads)
        else:
            FractalMatrix.path_to_fastaFiles(path, recursiv, filter_mikroSats=filter, threads=int(threads))
        if os.path.isfile(path):
            path = path + ".csv"
    for n in range(n_start, n_stop + 1):
        for k in range(k_start, n):
            if (k > k_stop):
                break
            FractalMatrix.path_to_markovPatternAnalyse(path, n, k, recursiv,log)

# cli()