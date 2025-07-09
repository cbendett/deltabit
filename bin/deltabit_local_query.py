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


parser = argparse.ArgumentParser(description='For all proteins in a query fasta, calculates their deltabits against an ingroup mtdb and outgroup mtdb. Outputs tsv.')
parser.add_argument('--query_faa', type=str, required=True, help='Path to the faa that will be the query for blasting')
parser.add_argument('--ingroup_mtdb', type=str, required=True, help='Path to the ingroup mtdb file')
parser.add_argument('--outgroup_mtdb', type=str, required=True, help='Path to the outgroup mtdb file')
parser.add_argument('--threads', type=int, default=1, help='How many threads to run the blast searches on')
parser.add_argument('--keep_intermediates', action='store_true', help='Add flag to not delete intermediate files made. Warning that these files can be quite large (into GB range)')
parser.add_argument('--output_name', type=str, required=True, help='Name for output files')
args = parser.parse_args()
print(args)


# In[ ]:


def ingroup_blast():
    #Define the command
    command = f'db2search -d {args.ingroup_mtdb} -q {args.query_faa} -a blastp --diamond -o ingroup_blast -c {args.threads} -m 3'
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
    command1 = f'db2search -d {args.outgroup_mtdb} -q {args.query_faa} -a blastp --diamond -o outgroup_blast -c {args.threads} -m 3'
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
    
    #Set qseqid as the index
    ingroup_df.set_index('qseqid', inplace=True)

    #Get ome of best hit
    ingroup_df['ingroup_ome'] = ingroup_df['sseqid'].str.split('_').str[0]

    #Sort by the bitscore (highest first)
    ingroup_df.sort_values(by=['bitscore'], ascending=False, inplace=True)

    #Keep only the highest bitscore entry for each qseqid
    ingroup_df = ingroup_df.groupby(ingroup_df.index).first()

    #Reset index
    ingroup_df = ingroup_df.reset_index(drop=False)

    #Keep only relevant info from ingroup_df (getting rid of extra blast info)
    ingroup_df_simple = ingroup_df[['qseqid', 'bitscore', 'ingroup_ome', 'sseqid']]
    ingroup_df_simple = ingroup_df_simple.rename(columns={
        'bitscore': 'ingroup_bitscore',
        'sseqid': 'ingroup_hit'
    })
    
    return(ingroup_df_simple)

ingroup_df = ingroup_consolidate('ingroup_blast.tsv')


# In[ ]:


def outgroup_consolidate(blast_table):
    #Read outgroup tsv into a dataframe and name columns
    outgroup_df = pd.read_csv(blast_table, sep='\t')
    outgroup_df.columns = ['qseqid', 'sseqid', 'pident', 'ppos', 'sstart', 'send', 'evalue', 'bitscore']

    #Set qseqid as the index
    outgroup_df.set_index('qseqid', inplace=True)

    #Get ome of best hit
    outgroup_df['outgroup_ome'] = outgroup_df['sseqid'].str.split('_').str[0]

    #Sort by the bitscore (highest first)
    outgroup_df.sort_values(by=['bitscore'], ascending=False, inplace=True)

    #Keep only the highest bitscore entry for each qseqid
    outgroup_df = outgroup_df.groupby(outgroup_df.index).first()

    #Reset index
    outgroup_df = outgroup_df.reset_index(drop=False)

    #Simplify outgroup_df for merging with ingroup_df
    outgroup_df_simple = outgroup_df[['qseqid', 'bitscore', 'sseqid', 'outgroup_ome']]
    outgroup_df_simple = outgroup_df_simple.rename(columns={
        'bitscore': 'outgroup_bitscore',
        'sseqid': 'outgroup_hit'
    })

    return(outgroup_df_simple)

outgroup_df = outgroup_consolidate('outgroup_blast.tsv')


# In[ ]:


#Merge ingroup_df and outgroup_df
merged_df = pd.merge(ingroup_df, outgroup_df, on='qseqid')

#Calculate delta_bit for each gene
merged_df['delta_bit'] = merged_df['ingroup_bitscore'] - merged_df['outgroup_bitscore']


# In[ ]:


merged_df.shape


# In[ ]:


#Export merged_df as tsv
merged_df.to_csv(f'{args.output_name}_deltabit.tsv', sep='\t', index=True)


# In[ ]:


# Create histogram of deltabits
sns.histplot(merged_df["delta_bit"], bins=30, kde=True)

# Label the axes
plt.xlabel("Deltabit")
plt.ylabel("Frequency")
plt.title("Histogram of Deltabit Values")
plt.tight_layout()

#Save the plot
plt.savefig(f'{args.output_name}_deltabit.png', format='png')

# Show the plot
plt.show()


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




