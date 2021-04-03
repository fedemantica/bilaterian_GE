#!/usr/bin/bash

#$1 reads the first argument from the command line

snakemake -j $1 --use-conda --cluster-config cluster.json --cluster "{cluster.qsub} -l h_vmem={cluster.ram},h_rt={cluster.time},disk={cluster.disk} -q {cluster.queue}" --keep-going --latency-wait 100
#--restart-times 2 --rerun-incomplete
