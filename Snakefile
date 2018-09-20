import os
import gzip
import itertools

# shortcuts
pj = os.path.join

# directories
data = config.get('data', 'data')
assemblies = pj(data, 'assemblies')
variants = pj(data, 'variants')
mutfunc = pj(data, 'mutfunc')
out = config.get('out', 'out')
ko = pj(out, 'hillenmeyer2008')
corr = pj(out, 'correlations')

# data files
raw = pj(data, 'ko_scores.tsv.gz')
rawref = pj(data, 'ko_scores_s288c.tsv.gz')
rawrep = pj(data, 'ko_scores_rep.tsv.gz')
rawsizes = pj(data, 'corrected.tsv.gz')
todrop = pj(data, 'S288C_to_drop.txt')
conditions = pj(data, 'conditions.tsv')
uniprot2gene = pj(data, 'uniprot2orf.tsv')
essential = pj(data, 'essentials.csv')
reactome_input = pj(data, 'reactomePathwaysScerevisiae.tsv')

# variables extracted from data file
strains = sorted({x.decode().rstrip().split('\t')[1]
                  for x in gzip.open(raw)} - {'strain'})

# output files
scores = pj(out, 'ko_scores.txt')
scoresref = pj(out, 'ko_scores_s288c.txt')
scoresrep = pj(out, 'ko_scores_rep.txt')
sizes = pj(out, 'sizes.txt')
fitness = pj(out, 'fitness.txt')
dups = pj(out, 'duplicates_correlation.tsv')
orth = pj(out, 'orthologs_correlation.tsv')
cond = pj(out, 'orthologs_conditions_correlation.tsv')
genes = pj(out, 'stratified_genes.tsv')
# conditions correlations
ccorrelations = pj(out, 'conditions_correlations.tsv')
# gene-gene correlations
scorrelations = [pj(corr, '%s.tsv' % x)
                 for x in strains]
pcorrelations = [pj(corr, '%s_%s.tsv' % (s1, s2))
                 for s1,s2 in itertools.combinations(strains, 2)]
# benchmarks
benchmarks = pj(out, 'benchmarks.tsv')
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
go_sets = pj(out, 'go.txt')
reactome = pj(out, 'reactome.txt')
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
# gene sets
gene_sets_tests = pj(out, 'gene_sets_tests.tsv')
# variants data
vcf = pj(variants, 'SGRP2-cerevisiae-freebayes-snps-Q30-GQ30.vcf.gz')
cvcf = pj(out, 'SGRP2-corrected.vcf.gz')
nvcf = pj(out, 'SGRP2-norm.vcf.gz')
vvcf = pj(out, 'SGRP2-variants.tsv')
pvcf = pj(out, 'SGRP2-plink.bed')
tvcf = pj(out, 'SGRP2.tsv')
# mutfunc
sift = pj(mutfunc, 'sift.tsv.gz')
foldx1 = pj(mutfunc, 'exp.tsv.gz')
foldx2 = pj(mutfunc, 'mod.tsv.gz')
mvcf = pj(out, 'SGRP2-mutfunc.tsv')

rule all:
  input:
    scores, scoresref, scoresrep, fitness,
    dups, scorrelations, pcorrelations,
    ccorrelations,
    benchmarks,
    orth, cond, genes,
    sdups, sorth, scond,
    kolog, koz, kopval,
    features, deviations,
    gaf, obo, goe,
    cpx, kegg, string,
    biogrid, biogrid_physical, biogrid_genetic,
    go_sets, reactome,
    mash, gene_sets_tests,
    mvcf, tvcf

rule fix_raw:
  input: raw, conditions, todrop
  output: scores
  shell: 'src/fix_raw {input} > {output}'

rule fix_rawref:
  input: rawref, conditions, todrop
  output: scoresref
  shell: 'src/fix_raw {input} > {output}'

rule fix_rawrep:
  input: rawrep, conditions, todrop
  output: scoresrep
  shell: 'src/fix_raw {input} > {output}'

rule fix_rawsizes:
  input: rawsizes, conditions, todrop
  output: sizes
  shell: 'src/fix_raw {input} > {output}'

rule fitness:
  input: sizes
  output: fitness
  shell: 'src/get_general_fitness {input} -r S288C -s Y55 YPS UWOP > {output}'

rule duplicates_correlation:
  input: scores
  output: dups
  shell: 'src/duplicates_correlation {input} > {output}'

rule conditions_correlation:
  input: scores
  output: ccorrelations
  shell: 'src/conditions_correlation {input} > {output}'

rule orthologs_correlation:
  input: scores
  output: orth
  shell: 'src/orthologs_correlation {input} > {output}'

rule orthologs_correlation_conditions:
  input: scores
  output: cond
  shell: 'src/orthologs_correlation_conditions {input} > {output}'

rule gene_correlations_single:
  input: scores
  output: scorrelations
  params: corr
  shell: 'src/get_genes_correlations {input} --out {params} --single'

rule gene_correlations:
  input: scores
  output: pcorrelations
  params: corr
  shell: 'src/get_genes_correlations {input} --out {params}'

rule benchmarking:
  input:
    a=scorrelations,
    b=cpx,
    c=kegg,
    d=biogrid_physical
  output: benchmarks
  params: corr
  shell: 'src/benchmark_correlations {params} {input.b} {input.c} {input.d} > {output}'

rule stratify_genes:
  input: scores
  output: genes
  shell: 'src/stratify_genes {input} > {output}'

rule duplicates_correlation_stratified:
  input:
    scores=scores,
    genes=genes
  output: pj('out', 'duplicates_correlation_{stratum}.tsv')
  shell: 'src/duplicates_correlation --genes {genes} --stratum {wildcards.stratum} {scores} > {output}'

rule orthologs_correaltion_stratified:
  input:
    scores=scores,
    genes=genes
  output: pj('out', 'orthologs_correlation_{stratum}.tsv')
  shell: 'src/orthologs_correlation --genes {genes} --stratum {wildcards.stratum} {scores} > {output}'

rule orthologs_correlation_conditions_stratified:
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
  input: gaf
  output: go_sets
  shell: 'cat {input} | src/gaf2txt > {output}'

rule:
  input: reactome_input
  output: reactome
  shell: 'cat {input} | src/reactome2txt > {output}'

rule gene_sets:
  input:
    scores,
    go_sets,
    cpx,
    reactome,
    kegg
  output:
    gene_sets_tests
  shell:
    'src/test_gene_sets {input} > {output}'

rule deviating_scores:
  input: scores, scoresrep
  output: deviations
  shell: 'src/get_deviating_scores {input} > {output}'

rule:
  input: scores, deviations 
  output: study, population
  shell: 'src/get_deviations_gene_sets {input} {output}'

rule deviating_genes_go_enrichment:
  input: obo, study, population, gaf
  output: goe
  shell: 'find_enrichment.py --obo {input} --method fdr_bh | grep "^GO" > {output}'

rule mash_sketches:
  input: assemblies
  output: mash_sketches
  shell:
    'mash sketch -s 10000 -o {output} {input}/*/assembly/genome.fa'

rule mash_distances:
  input: mash_sketches
  output: mash
  shell:
    'mash dist {input} {input} | src/square_mash > {output}'

rule download_vcf:
  output: vcf
  shell:
    'wget -O {output} "http://www.moseslab.csb.utoronto.ca/sgrp/data/SGRP2-cerevisiae-freebayes-snps-Q30-GQ30.vcf.gz"'

rule:
  input: vcf
  output: cvcf
  shell:
    'zcat {input} | src/correct_vcf | gzip > {output}'

rule:
  input: cvcf
  output: nvcf
  shell:  'bcftools norm -m - {input} | gzip > {output}'

rule:
  input: nvcf
  output: vvcf
  shell: 'src/annotate_vcf.R {input} {output}'

rule annotate_vcf:
  input:
    vcf=vvcf,
    sift=sift,
    exp=foldx1,
    mod=foldx2,
    conv=uniprot2gene
  output: mvcf
  shell: 'src/variants2mutfunc {input.vcf} {input.sift} {input.exp} {input.mod} --conversion {input.conv} > {output}'

rule:
  input: nvcf
  output: pvcf
  params: pj(out, os.path.split(pvcf)[1].split('.')[0])
  shell: 'plink --vcf {input} --out {params} --allow-extra-chr --make-bed'

rule matrix_vcf:
  input: pvcf
  output: tvcf
  params: pj(out, os.path.split(pvcf)[1].split('.')[0])
  shell: 'src/plink2df {params} > {output}' 
