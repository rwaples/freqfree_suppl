######   SETUP   ######

# set paths to the analysis programs,
# will need to be replaced your local installation
ANGSD="$HOME/programs/angsd/angsd"
realSFS="$HOME/programs/angsd/misc/realSFS"
IBS="$HOME/programs/angsd/misc/ibs"
SAMTOOLS="samtools"

# download the example data
wget http://popgen.dk/software/download/angsd/bams.tar.gz

# unzip/untar and index the bam files
tar xf bams.tar.gz
for i in bams/*.bam;do samtools index $i;done






######   realSFS METHOD   ######

# make a directory for the results
mkdir results_realsfs

# get the R script to parse the realSFS output
wget https://raw.githubusercontent.com/rwaples/freqfree_suppl/master/read_realSFS.R

# make a separate bam filelist for each indiviudal
# also create a SAMPLES array for use below
BAMS=./bams/*.bam
SAMPLES=()
for b in $BAMS; do
  # parse out the sample name
  base="$(basename -- $b)"
  sample="${base%%.mapped.*}"
  SAMPLES+=("$sample")
  echo $sample
  echo $b > ${sample}.filelist.ind
done

# For the realSFS method, one of the alleles at each site must be specified
# here we will use an ancestral state file
# download and index the ancestral state fasta file
wget http://popgen.dk/software/download/angsd/hg19ancNoChr.fa.gz
$SAMTOOLS faidx hg19ancNoChr.fa.gz


# run doSAF on each individual
for s in "${SAMPLES[@]}"; do
  $ANGSD -b ${s}.filelist.ind \
  -anc hg19ancNoChr.fa.gz \
  -minMapQ 30 -minQ 20 -GL 2 \
  -doSaf 1 -doDepth 1 -doCounts 1 \
  -out ${s}
done

# run realSFS on each pair of indiviudals
for i in {0..9}; do
  for j in {0..9}; do
    if (( i < j)); then
      sample1=${SAMPLES[i]}
      sample2=${SAMPLES[j]}
      $realSFS ${sample1}.saf.idx ${sample2}.saf.idx > ./results_realsfs/${sample1}_${sample2}.2dsfs
    fi
  done
done

# lets take a look at the results for a single pair of individuals
Rscript \
  -e "source('./read_realSFS.R')" \
  -e "res = read_realSFS('results_realsfs/smallNA06985_smallNA11830.2dsfs')" \
  -e "res['sample1'] = 'smallNA06985'; res['sample2'] = 'smallNA11830'" \
  -e "print(res[,c('sample1', 'sample2', 'nSites', 'Kin', 'R0', 'R1') ])"






######   IBS METHOD   ######

# make a directory for the results
mkdir results_IBS

# get the R script to parse the IBS output
wget https://raw.githubusercontent.com/rwaples/freqfree_suppl/master/read_IBS.R

## make a bam filelists containing all individuals
ls bams/*.bam > all.filelist

# make a genotype likelihood file (glf) containing all individuals
$ANGSD -b all.filelist \
  -minMapQ 30 -minQ 20 -GL 2 \
  -doGlf 1 \
  -out example

# run IBS, this will analyse each pair of indiviudals
$IBS -glf example.glf.gz \
  -model 0 \
  -nInd 10 -allpairs 1 \
  -outFileName results_IBS/ibs.model0.results

# look at the results from IBS
Rscript \
  -e "source('./read_IBS.R')" \
  -e "res = do_derived_stats(read_ibspair_model0('results_IBS/ibs.model0.results.ibspair'))" \
  -e "print(res[6,c('ind1', 'ind2', 'nSites', 'Kin', 'R0', 'R1') ])"

# the IBS method in ANGSD indexes individuals as they appear in the filelist
# (zero-indexed)
cat all.filelist
