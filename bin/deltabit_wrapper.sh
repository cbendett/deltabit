#!/bin/bash

### This is a wrapper bash script that coordinates the deltabit "deep" pipeline that performs deltabit on the entire contig of a specific captain gene, performs it on the buscos for a null distribution, identifies and extracts a region near the captain that has consecutive low deltabit scores and dubs it an htr, finds other genomes in the ingroup and outgroup that have clustered homologs of these htr genes, and produces a synteny plot to visualize it. 

set -e

############################ Step 0: Parse command line arguments ######################################

# Set default values for optional arguments
THREADS=16
CAPTAIN_FEATURE_TYPE="gene"
HOMOLOGS_TOP=10
HOMOLOGS_BITSCORE=40
HOMOLOGS_MINIMUM=3
HOMOLOGS_SORTED=20
CLUSTER_SIZE=20000

# Use getopt to parse long options
OPTS=$(getopt -o '' \
--long mycotools_db:,captain:,ingroup_mtdb:,path_to_deltabit:,captain_feature_type:,gff:,outgroup_mtdb_shallow:,outgroup_mtdb_deep:,threads:,tarred_busco_db:,faa:,busco_db_name:,prefix:,homologs_bitscore:,homologs_top:,homologs_minimum:,homologs_sorted:,cluster_size:,help \
-n 'deltabit_wrapper.sh' -- "$@")

if [ $? -ne 0 ]; then
  echo "Failed parsing options." >&2
  echo "Use --help to see usage." >&2
  exit 1
fi

eval set -- "$OPTS"

# Extract options
while true; do
  case "$1" in
    --help)
      echo "Usage: deltabit_wrapper.sh [OPTIONS]"
      echo
      echo "Required arguments:"
      echo "  --mycotools_db=PATH          Path to MycoTools database"
      echo "  --captain=PATH               Path to CAPTAIN"
      echo "  --ingroup_mtdb=PATH          Ingroup MycoTools database"
      echo "  --path_to_deltabit=PATH      Path to deltabit.py script"
      echo "  --gff=FILE                   GFF annotation file"
      echo "  --faa=FILE                   Protein FASTA file"
      echo "  --busco_db_name=NAME         Name of the BUSCO database"
      echo "  --prefix=STRING              Output prefix"
      echo "  --outgroup_mtdb_shallow=PATH Outgroup MycoTools database for whole contig and buscos"
      echo "  --outgroup_mtdb_deep=PATH    Outgroup MycoTools database for searching for other elements"
      echo "  --tarred_busco_db=FILE       Tarred BUSCO database"
      echo
      echo "Optional arguments:"
      echo "  --captain_feature_type=TYPE  Feature type in GFF (default: gene)"
      echo "  --threads=INT                Number of threads (default: 16)"
      echo "  --homologs_bitscore=INT      Homologs bitscore threshold (default: 40)"
      echo "  --homologs_top=INT           Top N homologs (default: 10)"
      echo "  --homologs_minimum=INT       Minimum homologs required (default: 3)"
      echo "  --homologs_sorted=INT        Sorted homologs to use (default: 20)"
      echo "  --cluster_size=INT           Cluster size (default: 20000)"
      echo "  --help                       Show this help message and exit"
      exit 0
      ;;
    --mycotools_db) MYCOTOOLS_DB="$2"; shift 2 ;;
    --captain) CAPTAIN="$2"; shift 2 ;;
    --ingroup_mtdb) INGROUP_MTDB="$2"; shift 2 ;;
    --path_to_deltabit) PATH_TO_DELTABIT="$2"; shift 2 ;;
    --captain_feature_type) CAPTAIN_FEATURE_TYPE="$2"; shift 2 ;;
    --gff) GFF="$2"; shift 2 ;;
    --outgroup_mtdb_shallow) OUTGROUP_MTDB_SHALLOW="$2"; shift 2 ;;
    --outgroup_mtdb_deep) OUTGROUP_MTDB_DEEP="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --tarred_busco_db) TARRED_BUSCO_DB="$2"; shift 2 ;;
    --faa) FAA="$2"; shift 2 ;;
    --busco_db_name) BUSCO_DB_NAME="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    --homologs_bitscore) HOMOLOGS_BITSCORE="$2"; shift 2 ;;
    --homologs_top) HOMOLOGS_TOP="$2"; shift 2 ;;
    --homologs_minimum) HOMOLOGS_MINIMUM="$2"; shift 2 ;;
    --homologs_sorted) HOMOLOGS_SORTED="$2"; shift 2 ;;
    --cluster_size) CLUSTER_SIZE="$2"; shift 2 ;;
    --) shift; break ;;
    *) echo "Invalid option: $1" >&2; echo "Use --help to see usage." >&2; exit 1 ;;
  esac
done

# Check for required arguments
missing_args=()

[ -z "$MYCOTOOLS_DB" ] && missing_args+=("--mycotools_db")
[ -z "$CAPTAIN" ] && missing_args+=("--captain")
[ -z "$INGROUP_MTDB" ] && missing_args+=("--ingroup_mtdb")
[ -z "$PATH_TO_DELTABIT" ] && missing_args+=("--path_to_deltabit")
[ -z "$GFF" ] && missing_args+=("--gff")
[ -z "$FAA" ] && missing_args+=("--faa")
[ -z "$BUSCO_DB_NAME" ] && missing_args+=("--busco_db_name")
[ -z "$PREFIX" ] && missing_args+=("--prefix")
[ -z "$OUTGROUP_MTDB_SHALLOW" ] && missing_args+=("--outgroup_mtdb_shallow")
[ -z "$OUTGROUP_MTDB_DEEP" ] && missing_args+=("--outgroup_mtdb_deep")
[ -z "$TARRED_BUSCO_DB" ] && missing_args+=("--tarred_busco_db")

if [ ${#missing_args[@]} -ne 0 ]; then
  echo "Error: Missing required arguments: ${missing_args[*]}" >&2
  echo "Use --help to see usage." >&2
  exit 1
fi

############## Step 1: Run deltabit.py on the entire contig of the specified captain gene ###################

# Make directory for Step 1
mkdir 01_whole_contig

# Set up mycotools 
mtdb -i "${MYCOTOOLS_DB}"

# Get ome from standard mycotools captain gene name
OME=$(echo "${CAPTAIN}" | awk -F'_' '{print $2}')

# Make ingroup mtdb that excludes genus
GENUS=$(mtdb $OME | awk '{print $2}')
grep -v "$GENUS" ${INGROUP_MTDB} > 01_whole_contig/ingroup.mtdb

# Run deltabit 
python "${PATH_TO_DELTABIT}"/deltabit.py --captain "${CAPTAIN}" --captain_feature_type "${CAPTAIN_FEATURE_TYPE}" --captain_db "${GFF}" --prox_dist 20000000000000000 --ingroup_mtdb 01_whole_contig/ingroup.mtdb --outgroup_mtdb "${OUTGROUP_MTDB_SHALLOW}" --threads "${THREADS}"

# Move results to 01_whole_contig
mv "${CAPTAIN}"_deltabit.png 01_whole_contig/"${CAPTAIN}"_whole_contig_deltabit.png
mv "${CAPTAIN}"_deltabit.tsv 01_whole_contig/"${CAPTAIN}"_whole_contig_deltabit.tsv

############# Step 2: Run deltabit_local_query.py on the buscos (and run busco itself) #######################

# Make directory for Step 2
mkdir 02_busco

# Set up busco dataset
mkdir 02_busco/lineages
tar -xzf "${TARRED_BUSCO_DB}" -C 02_busco/lineages/

# Run busco
busco -i "${FAA}" -l "${BUSCO_DB_NAME}" -o 02_busco/"${PREFIX}"_busco -m prot -c "${THREADS}" --download_path 02_busco/ --offline -f

# Get query faa from busco outputs
cat 02_busco/"${PREFIX}"_busco/run_"${BUSCO_DB_NAME}"/busco_sequences/single_copy_busco_sequences/*.faa > 02_busco/busco_query.faa

# Run deltabit with local query
python "${PATH_TO_DELTABIT}"/deltabit_local_query.py --query_faa 02_busco/busco_query.faa --ingroup_mtdb 01_whole_contig/ingroup.mtdb --outgroup_mtdb "${OUTGROUP_MTDB_SHALLOW}" --threads "${THREADS}" --output_name "${PREFIX}"_busco

# Move files to 02_busco
mv "${PREFIX}"_busco* 02_busco/

############## Step 3: Run deltabit_plot #######################################################################

# Make directory for Step 3 (and move into it)
mkdir 03_plot
cd 03_plot

# Run deltabit_plot.R
Rscript "${PATH_TO_DELTABIT}"/deltabit_plot.R -p "${PREFIX}" -d ../01_whole_contig/"${CAPTAIN}"_whole_contig_deltabit.tsv -b ../02_busco/"${PREFIX}"_busco_deltabit.tsv

# Move out of directory
cd ..

############## Step 4: BLAST sig_genes against the outgroup and find putative regions in other genomes #########

# Make directory for Step 4 (and move into it)
mkdir 04_blast_outgroup
cd 04_blast_outgroup 

# Make query fasta of sig_genes
acc2fa -i ../03_plot/"${PREFIX}"_sig_genes.txt > "${PREFIX}"_sig_genes.faa

# BLAST sig_genes against outgroup mtdb
db2search -a blastp --diamond -q "${PREFIX}"_sig_genes.faa -d "${OUTGROUP_MTDB_DEEP}" -c "${THREADS}" -m 3

# cat BLAST output reports
cat db*/reports/* > "${PREFIX}"_outgroup_blast.tsv

# Run deltabit_homologs.R
Rscript "${PATH_TO_DELTABIT}"/deltabit_homologs.R -p "${PREFIX}" -f "${PREFIX}"_outgroup_blast.tsv -b "${HOMOLOGS_BITSCORE}" -t "${HOMOLOGS_TOP}" -m "${HOMOLOGS_MINIMUM}" -s "${HOMOLOGS_SORTED}"

# Concatenate all homologs.txt files
cat *homologs.txt > "${PREFIX}"_homologs.txt

# Make gff of all homologs
acc2gff -i "${PREFIX}"_homologs.txt > "${PREFIX}"_homologs.gff

# Run deltabit_synteny.sh
"${PATH_TO_DELTABIT}"/deltabit_synteny.sh "${PREFIX}"_homologs.gff "${CLUSTER_SIZE}"

# Move out of directory
cd ..

############## Step 5: BLAST sig_genes against the ingroup and find putative regions in other genomes #########

# Make directory for Step 5 (and move into it)
mkdir 05_blast_ingroup
cd 05_blast_ingroup 

# Make query fasta of sig_genes
acc2fa -i ../03_plot/"${PREFIX}"_sig_genes.txt > "${PREFIX}"_sig_genes.faa

# BLAST sig_genes against ingroup mtdb
db2search -a blastp --diamond -q "${PREFIX}"_sig_genes.faa -d ../01_whole_contig/ingroup.mtdb -c "${THREADS}" -m 3

# cat BLAST output reports
cat db*/reports/* > "${PREFIX}"_ingroup_blast.tsv

# Run deltabit_homologs.R
Rscript "${PATH_TO_DELTABIT}"/deltabit_homologs.R -p "${PREFIX}" -f "${PREFIX}"_ingroup_blast.tsv -b "${HOMOLOGS_BITSCORE}" -t "${HOMOLOGS_TOP}" -m "${HOMOLOGS_MINIMUM}" -s "${HOMOLOGS_SORTED}"

# Concatenate all homologs.txt files
cat *homologs.txt > "${PREFIX}"_homologs.txt

# Make gff of all homologs
acc2gff -i "${PREFIX}"_homologs.txt > "${PREFIX}"_homologs.gff

# Run deltabit_synteny.sh
"${PATH_TO_DELTABIT}"/deltabit_synteny.sh "${PREFIX}"_homologs.gff "${CLUSTER_SIZE}"

# Move out of directory
cd ..

############# Step 6: Make synteny plot of HTR with other putative HTRs #######################################

# Make directory for Step 6 (and move into it)
mkdir 06_synteny
cd 06_synteny

# Make directory that will contain gbks
mkdir gbk

# Move gbks from 04_blast_outgroup and 05_blast_ingroup
cp ../04_blast_outgroup/*.gbk gbk/
cp ../05_blast_ingroup/*.gbk gbk/

# Make gbk of htr and move into gbk/
acc2gbk -i ../03_plot/"${PREFIX}"_htr.txt > "${PREFIX}".gbk
mv "${PREFIX}".gbk gbk/

# Count windows in gbk directory and save as variable for formatting
WINDOW_COUNT=$(ls gbk | wc -l)
SYNTENY_HEIGHT=$((WINDOW_COUNT * 90))

# Run deltabit_synteny.R
Rscript "${PATH_TO_DELTABIT}"/deltabit_synteny.R -f gbk -q "${PREFIX}" -p "${PREFIX}" -H "${SYNTENY_HEIGHT}"

