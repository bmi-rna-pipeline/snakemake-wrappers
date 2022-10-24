__author__ = "Seo Yoon Park"
__copyright__ = "Copyright 2022, Seo Yoon Park"
__email__ = "sypark02178@gmail.com"
# __license__ = "MIT"

import os

from snakemake.shell import shell

bam = snakemake.input.get("bam", "")
sam = snakemake.input.get("sam", "")
gtf = snakemake.input.get("gtf", "")

if bam:
    input_bam = " -f bam"
    input_string = bam
elif sam:
    input_bam = "-f sam"
    input_string = sam
else:
    raise Exception("Expected input.bam or input.sam, got neither.")


output_name = snakemake.output.countfiles

extra = snakemake.params.get("extra", "")
count = snakemake.params.get("count", "")

threads = snakemake.threads
log = snakemake.log_fmt_shell(stdout=True, stderr=True)
shell(
    "htseq-count -n {snakemake.threads} {extra} "
    "{input_bam} {input_string} "
    "--counts_output {output_name} "
    "> {log}"
)
