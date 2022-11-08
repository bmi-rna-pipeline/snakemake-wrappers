__author__ = "Seo Yoon Park"
__copyright__ = "Copyright 2022, Seo Yoon Park"
__email__ = ""
__license__ = ""


import tempfile
from pathlib import Path
from snakemake.shell import shell
from snakemake_wrapper_utils.samtools import get_samtools_opts

inputbam = snakemake.input.bamfile
outbam = snakemake.output.sortedout

samtools_opts = get_samtools_opts(snakemake)
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

awk = '{printf "%s", $0 ""; getline; print}'

with tempfile.TemporaryDirectory() as tmpdir:
    tmp_prefix = Path(tmpdir) / "samtools_fastq.sort_"

    shell(
        """cat <( samtools view -H {inputbam} )"""
        """<( samtools view -@ 12 {inputbam}"""
        """ | awk '{awk}'"""
        """ | sort -S 50G -T {tmp_prefix}"""
        """ | tr ' ' '\n' )"""
        """ | samtools view"""
        """ -@ 12 -bS -> {outbam}"""
        """ {log}"""
    )
