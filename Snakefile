import os
import itertools

# shortcuts
pj = os.path.join

# directories
data = config.get('data', 'data')
assemblies = pj(data, 'assemblies')
variants = pj(data, 'variants')
out = config.get('out', 'out')
ko = pj(out, 'hillenmeyer2008')
corr = pj(out, 'correlations')

# data files
raw = pj(data, 'ko_scores.tsv')
rawref = pj(data, 'ko_scores_s288c.tsv')
rawrep = pj(data, 'ko_scores_rep.tsv')
conditions = pj(data, 'conditions.tsv')

# variables extracted from data file
strains = sorted({x.rstrip().split('\t')[1]
                  for x in open(raw)} - {'strain'})

# output files
scores = pj(out, 'ko_scores.txt')
scoresref = pj(out, 'ko_scores_s288c.txt')
scoresrep = pj(out, 'ko_scores_rep.txt')
dups = pj(out, 'duplicates_correlation.tsv')
orth = pj(out, 'orthologs_correlation.tsv')
cond = pj(out, 'orthologs_conditions_correlation.tsv')
genes = pj(out, 'stratified_genes.tsv')
# gene-gene correaltions
scorrelations = [pj(corr, '%s.tsv' % x)
                 for x in strains]
pcorrelations = [pj(corr, '%s_%s.tsv' % (s1, s2))
                 for s1,s2 in itertools.combinations(strains, 2)]
#
# deprecated analysis
# COP/LLR scores and modules
cop = [pj(corr, '%s.cop.tsv' % x)
       for x in strains]
llr = [pj(corr, '%s.llr.tsv' % x)
       for x in strains]
llr_genes = pj(corr, 'unique_genes.txt')
modules = [pj(corr, '%s_%s.modules.tsv' % (s1, s2))
           for s1,s2 in itertools.combinations(strains, 2)
           if s1 == 'S288C'
           or s2 == 'S288C']
# go terms enrichment for conserved modules
mpopulation = pj(corr, 'population.txt')
gomodules = [pj(corr, '%s.go.tsv' % (s))
             for s in strains
             if s != 'S288C']
# end of deprecation
#
# deviating s-scores
deviations = pj(out, 'deviating.tsv')
# ko data (Hillenmeyer 2008)
kolog = pj(ko, 'lscores.tsv')
koz = pj(ko, 'zscores.tsv')
kopval = pj(ko, 'pvalues.tsv')
# SGD data
features = pj(out, 'SGD_features.tab')
gaf = pj(out, 'SGD_slim.tsv')
obo = pj(out, 'SGD_slim.obo')
# functional interactions
cpx = pj(out, 'complexes.cyc2008.txt')
kegg = pj(out, 'modules.kegg.txt')
string = pj(out, 'string.combined.800.txt')
biogrid = pj(out, 'biogrid.all.txt')
biogrid_physical = pj(out, 'biogrid.physical.txt')
biogrid_genetic = pj(out, 'biogrid.genetic.txt')
# variants data
snps = pj(variants, 'sgrp_Sc_SNPs.txt')
vcf = pj(variants, 'SGRP2-cerevisiae-freebayes-snps-Q30-GQ30.vcf.gz')
indels1 = pj(variants, 'indels', 'Sc_Ind_cr.txt')
indels2 = pj(variants, 'indels', 'Sc_Ind_ncr.txt')
# go terms enrichemnt
study = pj(out, 'deviating_study.txt')
population = pj(out, 'deviating_population.txt')
goe = pj(out, 'deviating_enrichemnt.tsv')
# stratified results
strata = ['g%d' % x for x in range(4)]
sdups = [pj('out', 'duplicates_correlation_%s.tsv' % x)
         for x in strata]
sorth = [pj('out', 'orthologs_correlation_%s.tsv' % x)
         for x in strata]
scond = [pj('out', 'orthologs_condition_correlation_%s.tsv' % x)
         for x in strata]
# mash distances
mash_sketches = pj(out, 'genomes.msh')
mash = pj(out, 'genome_distances.tsv')

rule all:
  input:
    scores, scoresref, scoresrep, dups,
    orth, cond, genes,
    sdups, sorth, scond,
    kolog, koz, kopval,
    features, deviations,
    gaf, obo, goe,
    cpx, kegg, string,
    biogrid, biogrid_physical, biogrid_genetic,
    mash, snps, vcf, indels1

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
  input:
    scores=scores,
    out=corr
  output: scorrelations
  shell: 'src/get_genes_correlations {input.scores} --out {input.out} --single'

rule:
  input:
    scores=scores,
    out=corr
  output: pcorrelations
  shell: 'src/get_genes_correlations {input.scores} --out {input.out}'

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

rule:
  output: gaf
  shell: 'curl --silent https://downloads.yeastgenome.org/curation/literature/go_slim_mapping.tab | src/slim2gaf > {output}'

rule:
  output: obo
  shell: 'wget -O {output} "http://www.geneontology.org/ontology/subsets/goslim_yeast.obo"'

rule:
  output: cpx
  shell: 'curl --silent "http://wodaklab.org/cyc2008/resources/CYC2008_complex.tab" | src/cyc2txt > {output}'

rule:
  output: kegg
  shell: 'src/get_kegg_modules sce | src/kegg2txt > {output}'

rule:
  output: string
  shell: 'wget --quiet --output-document - https://string-db.org/download/protein.links.detailed.v10.5/4932.protein.links.detailed.v10.5.txt.gz | zcat | src/string2txt --score 800 > {output}'

rule:
  output:
    biogrid=biogrid,
    biogrid_genetic=biogrid_genetic
  shell: 'wget -O tmp.zip https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-3.4.159/BIOGRID-ORGANISM-3.4.159.tab2.zip && unzip -j tmp.zip BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.4.159.tab2.txt -d . && cat BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.4.159.tab2.txt | src/biogrid2txt > {output.biogrid} && cat BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.4.159.tab2.txt | src/biogrid2txt --genetic > {output.biogrid_genetic} && rm tmp.zip && rm BIOGRID-ORGANISM-Saccharomyces_cerevisiae_S288c-3.4.159.tab2.txt'

rule:
  output: biogrid_physical
  shell: 'wget --quiet --output-document - https://downloads.thebiogrid.org/Download/BioGRID/Latest-Release/BIOGRID-MV-Physical-LATEST.tab2.zip | gunzip | src/biogrid2txt > {output}'

rule:
  input: scores, scoresrep
  output: deviations
  shell: 'src/get_deviating_scores {input} > {output}'

rule:
  input: scores, deviations 
  output: study, population
  shell: 'src/get_deviations_gene_sets {input} {output}'

rule:
  input: obo, study, population, gaf
  output: goe
  shell: 'find_enrichment.py --obo {input} --method fdr_bh | grep "^GO" > {output}'

rule:
  input:
    corr=pj(corr, 'S288C.tsv'),
    target=pj(corr, '{strain}.tsv'),
    cpx=biogrid_physical
  output: pj(corr, '{strain}.cop.tsv')
  shell:
    'src/get_cop_score {input.corr} {input.target} {input.cpx} --fraction 0.01 > {output}'

rule:
  input:
    corr=pj(corr, '{strain}.tsv'),
    cop=pj(corr, '{strain}.cop.tsv')
  output: pj(corr, '{strain}.llr.tsv')
  shell:
    'src/get_llr_score {input} > {output}'

rule:
  input: llr
  output: llr_genes
  shell: 'src/unique_interactions {input} > {output}'

rule:
  input:
    llr1=pj(corr, '{strain1}.llr.tsv'),
    llr2=pj(corr, '{strain2}.llr.tsv'),
    llr_genes=llr_genes
  output: pj(corr, '{strain1}_{strain2}.mergescore.tsv')
  shell: 'src/merge_score {input.llr1} {input.llr2} --subset {input.llr_genes} > {output}'

rule:
  input:
    scores=pj(corr, '{strain1}_{strain2}.mergescore.tsv'),
    llr1=pj(corr, '{strain1}.llr.tsv'),
    llr2=pj(corr, '{strain2}.llr.tsv')
  output: pj(corr, '{strain1}_{strain2}.modules.tsv')
  shell: 'src/merge_strains {input.scores} {input.llr1} {input.llr2} > {output}'

rule:
  input: llr_genes
  output: mpopulation
  shell: 'awk \'{{print $1"\\n"$2}}\' {input} | sort | uniq > {output}'

rule:
  input:
    modules=pj(corr, 'S288C_{strain}.modules.tsv'),
    obo=obo,
    pop=mpopulation,
    gaf=gaf
  output: pj(corr, '{strain}.go.tsv')
  shell: 'src/modules2go {input} {output} --minimum 6 --maximum 12 --spacer 100'

rule:
  input: assemblies
  output: mash_sketches
  shell:
    'mash sketch -s 10000 -o {output} {input}/*/assembly/genome.fa'

rule:
  input: mash_sketches
  output: mash
  shell:
    'mash dist {input} {input} | src/square_mash > {output}'

rule:
  output: snps
  shell:
    'wget -O {output} "http://www.moseslab.csb.utoronto.ca/sgrp/sgrp_Sc_SNPs.txt"'

rule:
  output: vcf
  shell:
    'wget -O {output} "http://www.moseslab.csb.utoronto.ca/sgrp/data/SGRP2-cerevisiae-freebayes-snps-Q30-GQ30.vcf.gz"'

rule:
  params: variants
  output: indels1
  shell:
    'wget -O {params}/indels.tar.gz "http://www.moseslab.csb.utoronto.ca/sgrp/sgrp_indels_apr2011.tar.gz" && cd {params} && tar -xvf indels.tar.gz && rm indels.tar.gz'
