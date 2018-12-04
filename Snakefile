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
ctodrop = pj(data, 'conditions_to_drop.txt')
conditions = pj(data, 'conditions.tsv')
rawnatural = pj(data, 'yeasts_natural.tsv.gz')
uniprot2gene = pj(data, 'uniprot2orf.tsv')
essential = pj(data, 'essentials.csv')
reactome_input = pj(data, 'reactomePathwaysScerevisiae.tsv')
genome = pj(data, 'reference.fasta')

# variables extracted from data file
strains = sorted({x.decode().rstrip().split('\t')[1]
                  for x in gzip.open(raw)} - {'strain'})

# output files
scores = pj(out, 'ko_scores.txt')
scoresref = pj(out, 'ko_scores_s288c.txt')
scoresrep = pj(out, 'ko_scores_rep.txt')
ascores = pj(out, 'ko_scores_annotated.txt')
wscores = pj(out, 'ko_scores_window.txt')
minsign = pj(out, 'ko_scores_minsignificance.txt')
sorted_conditions = pj(out, 'sorted_conditions.txt')
sorted_conditions_linkage = pj(out, 'sorted_conditions.linkage.gz')
sizes = pj(out, 'sizes.txt')
fitness = pj(out, 'fitness.txt')
natural = pj(out, 'yeasts_scores.txt')
natural_scores = pj(out, 'yeasts_sscores.txt')
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
s3correlation = pj(corr, 's3.tsv')
# benchmarks
benchmarks = pj(out, 'benchmarks.tsv')
s3benchmarks = pj(out, 'benchmarks_s3.tsv')
# deviating s-scores
deviations = pj(out, 'deviating.tsv')
# ko data (Hillenmeyer 2008)
kolog = pj(ko, 'lscores.tsv')
koz = pj(ko, 'zscores.tsv')
kopval = pj(ko, 'pvalues.tsv')
kolog1 = pj(ko, 'lscores_het.tsv')
koz1 = pj(ko, 'zscores_het.tsv')
kopval1 = pj(ko, 'pvalues_het.tsv')
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
strata = ['g%d' % x for x in range(3)]
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

# variants data - input and output
vcf = pj(variants, '1011Matrix.gvcf.gz')
rtab = pj(variants, 'genesMatrix_PresenceAbsence.tab.gz')
mat = pj(variants, '1011GWASMatrix.tar.gz')
matbed = pj(variants, '1011GWAS_matrix.bed')
matvcf = pj(variants, 'plink.vcf')
tree = pj(variants, '1011_matrix.tree.newick')
nvcf = pj(out, 'norm.vcf.gz')
similarity = pj(out, 'natural_similarity.tsv')
sickness = pj(data, 'sickness.tsv.gz')
avcf = pj(out, 'augmented.vcf.gz')
abed = pj(out, 'augmented.bed')
associations = pj(out, 'associations.tsv')
aassociations = pj(out, 'associations_annotated.tsv')
wassociations = pj(out, 'associations_window.tsv')
sgdbed = pj(out, 'SGD_features.bed')
sgdsortedbed = pj(out, 'SGD_sorted_features.bed')
genrichment = pj(out, 'gwas_enrichments.tsv')

rule all:
  input:
    scores, scoresref, scoresrep, fitness,
    natural, ascores, wscores, minsign,
    dups, scorrelations, s3correlation, pcorrelations,
    ccorrelations, sorted_conditions,
    sorted_conditions_linkage,
    benchmarks, s3benchmarks,
    orth, cond, genes,
    sdups, sorth, scond,
    kolog, koz, kopval,
    kolog1, koz1, kopval1,
    features, deviations,
    gaf, obo, goe,
    cpx, kegg, string,
    biogrid, biogrid_physical, biogrid_genetic,
    go_sets, reactome,
    mash, gene_sets_tests,
    aassociations, wassociations, genrichment

rule fix_raw:
  input: raw, conditions, todrop, ctodrop
  output: scores
  shell: 'src/fix_raw {input} > {output}'

rule annotate_scores:
  input: scores, sgdsortedbed
  output: ascores
  shell: 'src/annotate_scores {input} > {output}'

rule window_scores:
  input: genome, ascores
  output: wscores
  shell: 'src/ko_windows {input} --window 10000 > {output}'

rule fix_rawref:
  input: rawref, conditions, todrop, ctodrop
  output: scoresref
  shell: 'src/fix_raw {input} > {output}'

rule fix_rawrep:
  input: rawrep, conditions, todrop, ctodrop
  output: scoresrep
  shell: 'src/fix_raw {input} > {output}'

rule minimum_significance:
  input: scores
  output: minsign
  shell: 'src/get_minimum_significance {input} > {output}'

rule download_rawsizes:
  output: rawsizes
  shell: 'wget -O {output} "https://www.ebi.ac.uk/~marco/corrected.tsv.gz"'

rule fix_rawsizes:
  input: rawsizes, conditions, todrop, ctodrop
  output: sizes
  shell: 'src/fix_raw {input} > {output}'

rule sort_conditions:
  input: scores
  output:
    sc1=sorted_conditions,
    sc2=sorted_conditions_linkage
  params: 'S288C'
  shell: 'src/sort_conditions {input} --strain {params} --save-linkage {output.sc2} > {output.sc1}'

rule fitness:
  input: sizes
  output: fitness
  shell: 'src/get_general_fitness {input} -r S288C -s Y55 YPS UWOP > {output}'

rule fix_natural:
  input: rawnatural
  output: natural
  shell: 'src/fix_raw_natural {input} > {output}'

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

rule s3_correlations:
  input: scores
  output: s3correlation
  shell: 'src/get_genes_correlations_strains {input} > {output}'

rule benchmarking:
  input:
    a=scorrelations,
    b=cpx,
    c=kegg,
    d=biogrid_physical,
    e=minsign
  output: benchmarks
  params: corr
  shell:
    'src/benchmark_correlations {params} {input.b} {input.c} {input.d} --significance {input.e} --minimum-significance 1E-2 > {output}'

rule s3_benchmarking:
  input:
    a=s3correlation,
    b=cpx,
    c=kegg,
    d=biogrid_physical
  output: s3benchmarks
  params: corr
  shell: 'src/benchmark_correlations {params} {input.b} {input.c} {input.d} --strain s3 > {output}'

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
  output: kolog1
  shell: 'wget -O {output} "http://chemogenomics.stanford.edu/supplements/global/download/data/het.ratio_result_nm.pub"'

rule:
  output: koz1
  shell: 'wget -O {output} "http://chemogenomics.stanford.edu/supplements/global/download/data/het.z_result_nm.pub"'

rule:
  output: kopval1
  shell: 'wget -O {output} "http://chemogenomics.stanford.edu/supplements/global/download/data/het.z_tdist_pval_nm.pub"'

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
    'wget -O {output} "http://1002genomes.u-strasbg.fr/files/1011Matrix.gvcf.gz"'

rule download_rtab:
  output: rtab
  shell:
    'wget -O {output} "http://1002genomes.u-strasbg.fr/files/genesMatrix_PresenceAbsence.tab.gz"'

rule download_mat:
  output: mat
  shell:
    'wget -O {output} "http://1002genomes.u-strasbg.fr/files/1011GWASMatrix.tar.gz"'

rule unpack_mat:
  input: mat
  output: matbed
  params: variants
  shell: 'cd {params} && tar -xvf $(basename {input})'

rule normalize_vcf:
  input: vcf,
  output: nvcf
  shell:
    '''bcftools norm {input} -m - | sed 's/^chromosome//g' | bcftools view --min-ac 1 -q 0.05:minor | gzip > {output}'''

rule mat2vcf:
  input: matbed
  output: matvcf
  params: variants
  shell: 'cd {params} && plink --bfile $(basename {input} .bed) --recode vcf'

rule augment_vcf:
  input: nvcf, matvcf, rtab, sickness
  output: avcf
  shell: 'src/augment_vcf {input} | gzip > {output}'

rule vcf2bed:
  input: avcf
  output: abed
  params: out
  shell: 'plink --vcf {input} --maf 0.05 --recode12 --out {params}/$(basename {output} .bed) --allow-extra-chr --list-duplicate-vars suppress-first --make-bed'

rule similarity:
  input: tree
  output: similarity
  shell: 'python src/phylogeny_distance.py --lmm {input} > {output}'

rule get_sscores:
  input: natural
  output: natural_scores
  shell: 'src/get_sscores {input} > {output}'

rule associate:
  input: abed, natural_scores, similarity
  output: associations
  shell: 'src/limix_association {input} > {output}'

rule annotate_associations:
  input: associations, sgdsortedbed,
  output: aassociations
  shell: 'src/annotate_associations {input} > {output}'

rule window_associations:
  input: genome, ascores, aassociations
  output: wassociations
  shell: 'src/gwas_windows {input} --window 10000 > {output}'

rule sgd2bed:
  input: features
  output: sgdbed
  shell: 'src/sgd2bed {input} > {output}'

rule sort_bed:
  input: sgdbed
  output: sgdsortedbed
  shell: 'bedtools sort -i {input} > {output}'

rule gwas_enrichment:
  input: aassociations, genome, ascores, sgdsortedbed  
  output: genrichment
  shell:
    'src/gwas_enrichments {input} > {output}'
