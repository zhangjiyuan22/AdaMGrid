#!/bin/bash

mkdir -p logs_fullset

# amd-ep5-tmp:
# JOB_ID = 0 ... 1428
# total 1429 jobs
# at most 50 running simultaneously
sbatch \
    -p amd-ep5-tmp \
    --array=0-1428%50 \
    slurm.genmap_new_cluster_fullset

# amd-ep5:
# JOB_ID = 1429 ... 1999
# total 571 jobs
# at most 20 running simultaneously
sbatch \
    -p amd-ep5 \
    --array=1429-1999%20 \
    slurm.genmap_new_cluster_fullset