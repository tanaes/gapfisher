import click
import subprocess
from os.path import abspath

"""
This is the main Oecophylla script which handles launching of the entire
pipeline and installation of all necessary modules/environments.
Using a set of FASTA/FASTQ input files,
sample sheet (Illumina-specific) and tool parameters file,
it launches the Snakemake pipeline, either locally or on a supported
cluster system (Torque, Slurm).
Example usage:
==============
*NOT FUNCTIONAL AS OF NOW*
Installation:
-------------
oecophylla install
Execute Workflow:
-----------------
oecophylla workflow --input-dir ./inputs --sample-sheet sample.txt --params params.yaml --output-dir ./outputs
"""

@click.group()
def run():
    pass


def readfq(fp): # this is a generator function
    """
    FASTA/Q parser from Heng Li:
    https://github.com/lh3/readfq/blob/master/readfq.py
    """

    last = None # this is a buffer keeping the last unprocessed line
    while True: # mimic closure; is it a bad idea?
        if not last: # the first record or a record following a fastq
            for l in fp: # search for the start of the next record
                if l[0] in '>@': # fasta/q header line
                    last = l[:-1] # save this line
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp: # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+': # this is a fasta record
            yield name, ''.join(seqs), None # yield a fasta record
            if not last: break
        else: # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp: # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq): # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs); # yield a fastq record
                    break
            if last: # reach EOF before reading enough quality
                yield name, seq, None # yield a fasta record instead
                break

def clip_input(fasta, length):
    targets = {}

    for name, seq, _ in fasta:
        # if contig is less than 2x target region length,
        if len(seq) <= 2 * length:
            # just add it in 5' slot
            targets[name] = [seq,
                             (0, len(seq)),
                             None,
                             None]

        else:
            targets[name] = [seq[0:length],
                             (0, length),
                             seq[len(seq)-length:len(seq)],
                             (len(seq)-length, len(seq))]

    return(targets)

def format_clipped_fasta(targets):
    out = ''
    for target in targets:
        out += '>{0}_5\n{1}\n'.format(target,
                                      targets[target][0])
        if targets[target][2] is not None:
            out += '>{0}_3\n{1}\n'.format(target,
                                          targets[target][2])
    return(out)


def format_clipped_bed(targets):
    out = ''
    for target in targets:
        out += '{0}\t{1}\t{2}\n'.format(target,
                                        targets[target][1][0],
                                        targets[target][1][1])
        if targets[target][2] is not None:
            out += '{0}\t{1}\t{2}\n'.format(target,
                                            targets[target][3][0],
                                            targets[target][3][1])
    return(out)


def format_toml(fasta_in,
                targets,
                mmi,
                config,
                host,
                port):
    
    name=fasta_in
    control="false"
    min_chunks="0"
    max_chunks="inf"
    single_on="stop_receiving"
    multi_on="stop_receiving"
    single_off="unblock"
    multi_off="unblock"
    no_seq="proceed"
    no_map="unblock"

    targets_str = ''
    for target in targets:
        targets_str += '"{0}",{1},{2},-\n'.format(target,
                                                targets[target][1][0],
                                                targets[target][1][1])
        if targets[target][2] is not None:
            targets_str += '"{0}",{1},{2},+\n'.format(target,
                                                    targets[target][3][0],
                                                    targets[target][3][1])
        else:
            targets_str += '"{0}",{1},{2},+\n'.format(target,
                                                    targets[target][1][0],
                                                    targets[target][1][1])

    toml_out = TOML_BASE.format(config=config,
                                host=host,
                                port=port,
                                mmi_path=abspath(mmi),
                                name=name,
                                control=control,
                                min_chunks=min_chunks,
                                max_chunks=max_chunks,
                                targets=targets_str,
                                single_on=single_on,
                                multi_on=multi_on,
                                single_off=single_off,
                                multi_off=multi_off,
                                no_seq=no_seq,
                                no_map=no_map)
    return(toml_out)

def index_fasta(fasta_fp,
                mmi_fp):
    cmd = "minimap2 -x map-ont {0} -d {1}".format(fasta_fp,
                                                  mmi_fp)
    p1 = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    output = p1.communicate()[0]

    return(output)


@run.command()
# FASTA file(s) to index
@click.option('--input-fp', '-i', required=True, type=click.Path(exists=True),
               help='Input FASTA file with contigs')
# output: bed
@click.option('--bed', '-b', required=False, type=click.Path(),
              help='Output BED file with sequence targets')
# output: fasta
@click.option('--fasta', '-f', required=False, type=click.Path(),
              help='Output FASTA file of sequence targets')
# output: TOML
@click.option('--toml', '-t', required=False, type=click.Path(),
              help='Output TOML file specifying sequence targets')
# output: MMI
@click.option('--mmi', '-m', required=False, type=click.Path(),
              help='Output MiniMap Index of sequence target')
# param: length
@click.option('--length', '-l', required=False, type=click.INT, default=2000,
              help='Length of target region around contig ends')
# param: caller config
@click.option('--config', '-c', required=False, type=click.STRING,
              default='dna_r9.4.1_450bps_hac',
              help='ONT basecalling config to use')
# param: caller host
@click.option('--host-ip', required=False, type=click.STRING,
              default='127.0.0.1',
              help='IP address to basecaller host')
# param: caller port
@click.option('--host-port', required=False, type=click.STRING,
              default='5555',
              help='Host basecaller port')


def winnow(input_fp, bed, fasta, toml, mmi, length, config, host_ip, host_port):

    # check inputs
    if all(v is None for v in (bed, fasta, toml)):
        raise ValueError("Must select at least one of "
                         "--bed, --fasta, or --toml")

    if toml is not None and mmi is None:
        raise ValueError("If selecting --toml, must also"
                         "specify --mmi")

    # read input fasta
    with open(input_fp, 'r') as f:
        fasta_in = readfq(f)

        # iterate over contigs
        targets = clip_input(fasta_in, length)

    # if output fasta:
    if fasta is not None:
        with open(fasta, 'w') as f_o:
            fasta_out = format_clipped_fasta(targets)
            f_o.write(fasta_out)


    # if output bed:
    if bed is not None:
        with open(bed, 'w') as b_o:
            bed_out = format_clipped_bed(targets)
            b_o.write(bed_out)

    # if output TOML:
    if toml is not None:
        with open(toml, 'w') as t_o:
            toml_out = format_toml(input_fp,
                                   targets,
                                   mmi,
                                   config,
                                   host_ip,
                                   host_port)
            t_o.write(toml_out)

    # if index fasta input with MiniMap2:
    if mmi is not None:
        log = index_fasta(input_fp, mmi)

# - 

# Utilities:
# - simulate results -- on graph?
# - use bandage + plot graph
# - create BED 
# - 




# read in fasta file
# identify read ends
# create bed


TOML_BASE = (
            '[caller_settings]\n'
            'config_name = "{config}"\n'
            'host = "{host}"\n'
            'port = "{port}"\n'
            '\n'
            '[conditions]\n'
            'reference = "{mmi_path}"\n'
            '\n'
            '[conditions.0]\n'
            'name = "{name}"\n'
            'control = {control}\n'
            'min_chunks = {min_chunks}\n'
            'max_chunks = {max_chunks}\n'
            'targets = [{targets}]\n'
            'single_on = "{single_on}"\n'
            'single_off = "{single_off}"\n'
            'multi_on = "{multi_on}"\n'
            'multi_off = "{multi_off}"\n'
            'no_seq = "{no_seq}"\n'
            'no_map = "{no_map}"\n'
            )