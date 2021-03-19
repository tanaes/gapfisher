import pytest
import logging
from io import StringIO
from unittest import TestCase
from tempfile import TemporaryDirectory
from os.path import join, isfile
from gapfisher.gapfisher import (readfq, clip_input, format_clipped_fasta,
                                 format_clipped_bed, format_toml, index_fasta)

def test_readfq():
    fasta_f = StringIO(FASTA_IN)
    fasta = readfq(fasta_f)

    n, s, _ = next(fasta)

    assert(n == 'contig1')
    assert(s == 'AGCCCTTAGTCGCGCTACGATCGATCGATC'
                'AACAGTCGACGATCGTCGATCGCGACGATC')

def test_clip_input():
    fasta_f = StringIO(FASTA_IN)
    fasta = readfq(fasta_f)
    targets_obs = clip_input(fasta, 60)

    TestCase().assertDictEqual(targets_exp, targets_obs)

def test_format_clipped_fasta():
    obs_fasta_out = format_clipped_fasta(targets_exp)
    assert(obs_fasta_out == FASTA_OUT)

def test_format_clipped_bed():
    obs_bed_out = format_clipped_bed(targets_exp)
    assert(obs_bed_out == BED_OUT)

def test_format_toml():
    obs_toml_out = format_toml('fasta_fp.fasta',
                               targets_exp,
                               '/path/to/reference.mmi',
                               'dna_r9.4.1_450bps_hac',
                               '127.0.0.1',
                               '5555')

    assert(obs_toml_out == TOML_OUT)

def test_index_fasta():
    with TemporaryDirectory() as tempdir:
        fasta_fp = join(tempdir, 'fasta.fa')
        with open(fasta_fp, 'w') as f:
            f.write(FASTA_IN)
        mmi_fp = join(tempdir, 'fasta.mmi')

        out = index_fasta(fasta_fp, mmi_fp)

        assert(isfile(mmi_fp))



FASTA_IN = """>contig1
AGCCCTTAGTCGCGCTACGATCGATCGATC
AACAGTCGACGATCGTCGATCGCGACGATC

>contig2
AGCCCTTAGTCGCGCTACGATCGATCGATC
AACAGTCGACGATCGTCGATCGCGACGATC
AGCCCTTAGTCGCGCTACGATCGATCGATC
AACAGTCGACGATCGTCGATCGCGACGATC
AGCCCTTAGTCGCGCTACGATCGATCGATC
AACAGTCGACGATCGTCGATCGCGACGATC

>contig3 foo
AGCCCTTAGTCGCGCTACGATCGATCGATC
AACAGTCGACGATCGTCGATCGCGACGATC
AGCCCTTAGTCGCGCTACGATCGATCGATC
"""

FASTA_OUT = """>contig1_5
AGCCCTTAGTCGCGCTACGATCGATCGATCAACAGTCGACGATCGTCGATCGCGACGATC
>contig2_5
AGCCCTTAGTCGCGCTACGATCGATCGATCAACAGTCGACGATCGTCGATCGCGACGATC
>contig2_3
AGCCCTTAGTCGCGCTACGATCGATCGATCAACAGTCGACGATCGTCGATCGCGACGATC
>contig3_5
AGCCCTTAGTCGCGCTACGATCGATCGATCAACAGTCGACGATCGTCGATCGCGACGATCAGCCCTTAGTCGCGCTACGATCGATCGATC
"""

BED_OUT = """contig1\t0\t60
contig2\t0\t60
contig2\t120\t180
contig3\t0\t90
"""

TOML_OUT = """[caller_settings]
config_name = "dna_r9.4.1_450bps_hac"
host = "127.0.0.1"
port = "5555"

[conditions]
reference = "/path/to/reference.mmi"

[conditions.0]
name = "fasta_fp.fasta"
control = false
min_chunks = 0
max_chunks = inf
targets = ["contig1",0,60,-
"contig1",0,60,+
"contig2",0,60,-
"contig2",120,180,+
"contig3",0,90,-
"contig3",0,90,+
]
single_on = "stop_receiving"
single_off = "unblock"
multi_on = "stop_receiving"
multi_off = "unblock"
no_seq = "proceed"
no_map = "unblock"
"""

targets_exp = {'contig1': ['AGCCCTTAGTCGCGCTACGATCGATCGATC'
                               'AACAGTCGACGATCGTCGATCGCGACGATC',
                               (0, 60),
                               None,
                               None],
                   'contig2': ['AGCCCTTAGTCGCGCTACGATCGATCGATC'
                               'AACAGTCGACGATCGTCGATCGCGACGATC',
                               (0, 60),
                               'AGCCCTTAGTCGCGCTACGATCGATCGATC'
                               'AACAGTCGACGATCGTCGATCGCGACGATC',
                               (120, 180)],
                   'contig3': ['AGCCCTTAGTCGCGCTACGATCGATCGATC'
                               'AACAGTCGACGATCGTCGATCGCGACGATC'
                               'AGCCCTTAGTCGCGCTACGATCGATCGATC',
                               (0, 90),
                               None,
                               None]}