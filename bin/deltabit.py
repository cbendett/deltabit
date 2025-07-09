#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#For a given gene, identify the proximal genes, save their distance from the original gene. 
#Blast each proximal gene against an ingroup database and an outgroup database. 
#Use blast results to calculate alien index for each gene. 
#Output table that has gene, alien index, and distance from ‘captain.’


# In[ ]:


import argparse
import pandas as pd
import gffpandas.gffpandas as gffpd
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import subprocess
import seaborn as sns
import matplotlib.pyplot as plt
import os
import shutil


# In[ ]:


parser = argparse.ArgumentParser(description='For a given gene (captain), identifies proximal genes and calculates their alien index against designated ingroup and outgroup datasets. Outputs a TSV of all the proximal genes with their names, alien indices, and distances from the captain')
parser.add_argument('--captain', type=str, required=True, help='ID of the captain gene. Must match ID in database gff3 used for finding proximal genes. Do not refer to it by something else like name or alias.')
parser.add_argument('--captain_feature_type', type=str, default='gene', help='The type of feature the captain is according to your gff3. Default is gene.')
parser.add_argument('--captain_db', type=str, required=True, help='Path to the assembly gff3 that includes the captain and will be used to find proximal genes')
parser.add_argument('--prox_dist', type=int, default=100, help='How far from the captain to search for proximal genes in KB. Default is 100. If captain is toward end or beginning of scaffold, might not be able to search the entire distance.')
parser.add_argument('--ingroup_mtdb', type=str, required=True, help='Path to the ingroup mtdb file')
parser.add_argument('--outgroup_mtdb', type=str, required=True, help='Path to the outgroup mtdb file')
parser.add_argument('--threads', type=int, default=1, help='How many threads to run the blast searches on')
parser.add_argument('--qc_threshold', type=float, default=0, help='Will filter out any blast hits with a query cover less than the specified value (from 0 to 1). Default is 0.') 
parser.add_argument('--keep_intermediates', action='store_true', help='Add flag to not delete intermediate files made. Warning that these files can be quite large (into GB range)')
args = parser.parse_args()


# In[ ]:


#Define assembly gff path
assembly_gff_file = f'{args.captain_db}'
assembly_gff_df = gffpd.read_gff3(assembly_gff_file)
gene_row = assembly_gff_df.get_feature_by_attribute('ID', [f'{args.captain}']).filter_feature_of_type([f'{args.captain_feature_type}'])
captain_scaffold = gene_row.df['seq_id'].iloc[0]
captain_pos = (gene_row.df['start'].iloc[0] + gene_row.df['end'].iloc[0]) / 2


# In[ ]:


def extract_genes(gff_file, scaffold, captain_pos, search_range):
    #Load the GFF file
    gff = gffpd.read_gff3(gff_file)
    df = gff.df

     # Check if required columns exist in the DataFrame
    required_columns = ['seq_id', 'type', 'start', 'end', 'attributes']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns in the GFF file: {', '.join(missing_columns)}")

    #Filter by scaffold, type (gene), and base pair range
    genes_in_range = df[
        (df['seq_id'] == scaffold) &
        (df['type'] == 'gene') &
        (df['start'] >= (captain_pos - search_range)) &
        (df['end'] <= (captain_pos + search_range))
    ].copy()  #Copy to avoid SettingWithCopy warnings

    #Calculate distance from captain for each gene and add it as a new column
    genes_in_range['gene_midpoint'] = (genes_in_range['start'] + genes_in_range['end']) / 2
    genes_in_range['distance_from_captain'] = genes_in_range['gene_midpoint'] - captain_pos

    #Calculate gene length and add it as a new column
    genes_in_range['gene_length'] = genes_in_range['end'] - genes_in_range['start']

    #Extract gene name and make new column
    genes_in_range['gene_name'] = genes_in_range['attributes'].str.extract(r'Alias=([^;]+)')

    #Write gene names to new text file
    gene_names = genes_in_range['gene_name'].dropna()
    with open('proximal_genes.txt', 'w') as f:
        f.write('\n'.join(gene_names))

    return genes_in_range

genes_in_range = extract_genes(assembly_gff_file, captain_scaffold, captain_pos, args.prox_dist)


# In[ ]:


def acc2fa():
    #Define the command
    command = 'acc2fa -i proximal_genes.txt > proximal_genes.fa'

    #Run the command using subprocess
    try:
        subprocess.run(command, shell=True, check=True)
        print('Successful conversion from list to fasta.')
    except subprocess.CalledProcessError as e:
        print(f'An error in converting list to fasta occurred: {e}')

acc2fa()


# In[ ]:


def ingroup_blast():
    #Define the command
    command = f'db2search -d {args.ingroup_mtdb} -q proximal_genes.fa -a blastp --diamond -o ingroup_blast -c {args.threads}'
    command2 = 'cat ingroup_blast/reports/*.tsv > ingroup_blast.tsv'

    #Run command1 using subprocess
    try:
        subprocess.run(command, shell=True, check=True)
        print('Successful blast of proximal genes against ingroup mtdb')
    except subprocess.CalledProcessError as e:
        print(f'An error occurred in blasting the proximal genes against the ingroup mtdb: {e}')

    #Run command2 using subprocess
    try:
        subprocess.run(command2, shell=True, check=True)
        print('Successful concatenation of ingroup blast tsvs')
    except subprocess.CalledProcessError as e:
        print(f'An error occurred in concatenating the ingroup blast tsvs: {e}')

ingroup_blast()


# In[ ]:


def outgroup_blast():
    #Define the command
    command1 = f'db2search -d {args.outgroup_mtdb} -q proximal_genes.fa -a blastp --diamond -o outgroup_blast -c {args.threads}'
    command2 = 'cat outgroup_blast/reports/*.tsv > outgroup_blast.tsv'

    #Run command1 using subprocess
    try:
        subprocess.run(command1, shell=True, check=True)
        print('Successful blast of proximal genes against outgroup mtdb')
    except subprocess.CalledProcessError as e:
        print(f'An error occurred in blasting the proximal genes against the outgroup mtdb: {e}')

    #Run command2 using subprocess
    try:
        subprocess.run(command2, shell=True, check=True)
        print('Successful concatenation of outgroup blast tsvs')
    except subprocess.CalledProcessError as e:
        print(f'An error occurred in concatenating the outgroup blast tsvs: {e}')

outgroup_blast()


# In[ ]:


def ingroup_consolidate(blast_table):
    #Read ingroup tsv into a dataframe and name columns
    ingroup_df = pd.read_csv(blast_table, sep='\t')
    ingroup_df.columns = ['qseqid', 'sseqid', 'pident', 'ppos', 'sstart', 'send', 'evalue', 'bitscore']

    #Merge ingroup_df with genes_in_range
    ingroup_df = pd.merge(ingroup_df, genes_in_range, left_on='qseqid', right_on='gene_name', how='inner')

    #Add a new column of the hit length
    ingroup_df['hit_length'] = ingroup_df['send'] - ingroup_df['sstart']

    #Add a new column of the query cover
    ingroup_df['query_cover'] = ingroup_df['hit_length'] / ingroup_df['gene_length']

    #Filter out hits with low query cover
    ingroup_df = ingroup_df[ingroup_df['query_cover'] >= args.qc_threshold]
    
    #Set qseqid as the index
    ingroup_df.set_index('qseqid', inplace=True)

    #Get ome of best hit
    ingroup_df['ingroup_ome'] = ingroup_df['sseqid'].str.split('_').str[0]

    #Sort by the bitscore (highest first)
    ingroup_df.sort_values(by=['bitscore'], ascending=False, inplace=True)

    #Keep only the highest bitscore entry for each qseqid
    ingroup_df = ingroup_df.groupby(ingroup_df.index).first()

    #Keep only relevant info from ingroup_df (getting rid of extra blast info)
    ingroup_df_simple = ingroup_df[['gene_name', 'seq_id', 'start', 'end', 'gene_midpoint', 'distance_from_captain', 'gene_length', 'attributes', 'bitscore', 'query_cover', 'ingroup_ome', 'sseqid']].reset_index(drop=True)
    ingroup_df_simple = ingroup_df_simple.rename(columns={
        'bitscore': 'ingroup_bitscore',
        'query_cover': 'ingroup_query_cover',
        'sseqid': 'ingroup_hit'
    })
    
    return(ingroup_df_simple)

ingroup_df = ingroup_consolidate('ingroup_blast.tsv')


# In[ ]:


def outgroup_consolidate(blast_table):
    #Read outgroup tsv into a dataframe and name columns
    outgroup_df = pd.read_csv(blast_table, sep='\t')
    outgroup_df.columns = ['qseqid', 'sseqid', 'pident', 'ppos', 'sstart', 'send', 'evalue', 'bitscore']

    #Merge outgroup_df with genes_in_range
    outgroup_df = pd.merge(outgroup_df, genes_in_range, left_on='qseqid', right_on='gene_name', how='inner')

    #Add a new column of the hit length
    outgroup_df['hit_length'] = outgroup_df['send'] - outgroup_df['sstart']

    #Add a new column of the query cover
    outgroup_df['query_cover'] = outgroup_df['hit_length'] / outgroup_df['gene_length']

    #Filter out hits with low query cover
    outgroup_df = outgroup_df[outgroup_df['query_cover'] >= args.qc_threshold]
    
    #Set qseqid as the index
    outgroup_df.set_index('qseqid', inplace=True)

    #Get ome of best hit
    outgroup_df['outgroup_ome'] = outgroup_df['sseqid'].str.split('_').str[0]

    #Sort by the bitscore (highest first)
    outgroup_df.sort_values(by=['bitscore'], ascending=False, inplace=True)

    #Keep only the highest bitscore entry for each qseqid
    outgroup_df = outgroup_df.groupby(outgroup_df.index).first()

    #Simplify outgroup_df for merging with ingroup_df
    outgroup_df_simple = outgroup_df[['gene_name', 'bitscore', 'query_cover', 'sseqid', 'outgroup_ome']]
    outgroup_df_simple = outgroup_df_simple.rename(columns={
        'bitscore': 'outgroup_bitscore',
        'query_cover': 'outgroup_query_cover',
        'sseqid': 'outgroup_hit'
    })

    return(outgroup_df_simple)

outgroup_df = outgroup_consolidate('outgroup_blast.tsv')


# In[ ]:


#Merge ingroup_df and outgroup_df
merged_df = pd.merge(ingroup_df, outgroup_df, on='gene_name')

#Calculate delta_bit for each gene
merged_df['delta_bit'] = merged_df['ingroup_bitscore'] - merged_df['outgroup_bitscore']


# In[ ]:


#Export merged_df as tsv
merged_df.to_csv(f'{args.captain}_deltabit.tsv', sep='\t', index=True)


# In[ ]:


#Create basic plot of distance from captain vs delta_bit

#Create scatterplot
plt.figure(figsize=(10, 6))
sns.scatterplot(data=merged_df, x='distance_from_captain', y='delta_bit')

# Set y-axis to symmetric log scale
plt.yscale('symlog')

# Add labels and title
plt.xlabel('Distance from captain (bp)')
plt.ylabel('delta bit')
plt.title(f'{args.captain} delta bit plot')
plt.tight_layout()

#Save the plot
plt.savefig(f'{args.captain}_deltabit.png', format='png')


# In[ ]:


#Delete large intermediate files

#Paths to intermediates
files_to_delete = ['proximal_genes.fa', 'ingroup_blast.tsv', 'outgroup_blast.tsv']
folders_to_delete = ['ingroup_blast', 'outgroup_blast']

#Check for --keep_intermediates flag
if not args.keep_intermediates:
    #Delete files
    for file in files_to_delete:
        try:
            os.remove(file)
            print(f'Deleted: {file}')
        except FileNotFoundError:
            print(f'File not found: {file}')
        except Exception as e:
            print(f'Error deleting {file}: {e}')

    #Delete folders
    for folder in folders_to_delete:
        try:
            shutil.rmtree(folder)
            print(f'Deleted folder: {folder}')
        except FileNotFoundError:
            print(f'Folder not found: {folder}')
        except Exception as e:
            print(f'Error deleting folder {folder}: {e}')


# In[ ]: