#!/bin/bash

python3 ../../src/correct_raw_sizes.py sizes.txt.gz emap_info/ | gzip > corrected.tsv.gz
