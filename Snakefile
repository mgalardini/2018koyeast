import os

# shortcuts
pj = os.path.join

# directories
data = config.get('data', 'data')
out = config.get('out', 'out')
ko = pj(out, 'hillenmeyer2008')

# data files
raw = pj(data, 'ko_scores.tsv')
rawref = pj(data, 'ko_scores_s288c.tsv')
rawrep = pj(data, 'ko_scores_rep.tsv')
conditions = pj(data, 'conditions.tsv')

# output files
scores = pj(out, 'ko_scores.txt')
scoresref = pj(out, 'ko_scores_s288c.txt')
scoresrep = pj(out, 'ko_scores_rep.txt')
dups = pj(out, 'duplicates_correlation.tsv')
orth = pj(out, 'orthologs_correlation.tsv')
cond = pj(out, 'orthologs_conditions_correlation.tsv')
genes = pj(out, 'stratified_genes.tsv')
# ko data (Hillenmeyer 2008)
kolog = pj(ko, 'lscores.tsv')
koz = pj(ko, 'zscores.tsv')
kopval = pj(ko, 'pvalues.tsv')
# SGD data
features = pj(out, 'SGD_features.tab')

# stratified results
strata = ['g%d' % x for x in range(4)]
sdups = [pj('out', 'duplicates_correlation_%s.tsv' % x)
         for x in strata]
sorth = [pj('out', 'orthologs_correlation_%s.tsv' % x)
         for x in strata]
scond = [pj('out', 'orthologs_condition_correlation_%s.tsv' % x)
         for x in strata]

rule all:
  input:
    scores, scoresref, scoresrep, dups,
    orth, cond, genes,
    sdups, sorth, scond,
    kolog, koz, kopval,
    features

rule:
  input: raw, conditions
  output: scores
  shell: 'src/fix_conditions {input} > {output}'

rule:
  input: rawref, conditions
  output: scoresref
  shell: 'src/fix_conditions {input} > {output}'

rule:
  input: rawrep, conditions
  output: scoresrep
  shell: 'src/fix_conditions {input} > {output}'

rule:
  input: scores
  output: dups
  shell: 'src/duplicates_correlation {input} > {output}'

rule:
  input: scores
  output: orth
  shell: 'src/orthologs_correlation {input} > {output}'

rule:
  input: scores
  output: cond
  shell: 'src/orthologs_correlation_conditions {input} > {output}'

rule:
  input: scores
  output: genes
  shell: 'src/stratify_genes {input} > {output}'

rule:
  input:
    scores=scores,
    genes=genes
  output: pj('out', 'duplicates_correlation_{stratum}.tsv')
  shell: 'src/duplicates_correlation --genes {genes} --stratum {wildcards.stratum} {scores} > {output}'

rule:
  input:
    scores=scores,
    genes=genes
  output: pj('out', 'orthologs_correlation_{stratum}.tsv')
  shell: 'src/orthologs_correlation --genes {genes} --stratum {wildcards.stratum} {scores} > {output}'

rule:
  input:
    scores=scores,
    genes=genes
  output: pj('out', 'orthologs_condition_correlation_{stratum}.tsv')
  shell: 'src/orthologs_correlation_conditions --genes {genes} --stratum {wildcards.stratum} {scores} > {output}'

rule:
  output: kolog
  shell: 'wget -O {output} "http://chemogenomics.stanford.edu/supplements/global/download/data/hom.ratio_result_nm.pub"'

rule:
  output: koz
  shell: 'wget -O {output} "http://chemogenomics.stanford.edu/supplements/global/download/data/hom.z_result_nm.pub"'

rule:
  output: kopval
  shell: 'wget -O {output} "http://chemogenomics.stanford.edu/supplements/global/download/data/hom.z_tdist_pval_nm.pub"'

rule:
  output: features
  shell: 'wget -O {output} "https://downloads.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab"'