#!/bin/bash

cp ~/leia/2017yeast/20170928/emap/all_genes/ko_scores.tsv data/ko_scores.tsv
wget -O data/conditions.tsv "https://docs.google.com/spreadsheets/d/14i_5lWftFKMTwS4dOWULescrMCizwndyP-k6Y-P_SoI/export?gid=54998998&format=tsv"
