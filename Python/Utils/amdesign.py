""" Module for tools to design Antigen Arrays - currently mainly to generate overlapping peptide sets from a given set of antigens"""
import pandas as pd
import matplotlib.mlab as mlab
import numpy as np
import pylab
import itertools
from Bio import SeqIO
import ipdb as ipdb

# This is a test to see if commit works.

__all__ = ['create_am_peptide_set_from_fasta']


def create_am_peptide_set_from_fasta(filename, pep_length=20, pep_overlap=15, 
                                    insertion_positions=None, sort_method=None, output_filename=None):
    """create a peptide set from the given set of fasta files for designing peptide anitgen microarrays"""

    if output_filename is None:
        output_filename = "".join([filename[0:filename.rfind('.')], '_peptides_L', str(pep_length), '_O', str(pep_overlap), '.txt'])  

    pep_df = None
    columns = ['Name', 'begInd', 'endInd', 'Sequence']
    pep_data = {}
    fasta_sequences = SeqIO.parse(open(filename),'fasta')
    for fasta in fasta_sequences:
        name, sequence = fasta.id, fasta.seq.tostring()
        seq_len = len(sequence) 
        # handle insertions to correct for overlap. For each insertion start the beginning indices one step backwards to
        # allow for matches on all other positions to be retained (prevents large unrequired peptide sets due to insertions)

        pep_data['begInd'] = np.arange(0, len(sequence), pep_length-pep_overlap)
        pep_data['endInd'] = pep_data['begInd'] + pep_length
        pep_data['endInd'][pep_data['endInd'] >= seq_len] = seq_len-1
        pep_data['pep_names'] = []
        peptides = []

        if insertion_positions is not None:
            for ind in insertion_positions:
                pep_data['begInd'][np.where(pep_dat['begInd'] > ind)] += 1
    
        for i, inds in enumerate(zip(pep_data['begInd'], pep_data['endInd'])):
            peptides.append(sequence[inds[0]:inds[1]])
            pep_data['pep_names'].append("".join([name[0:name.rfind('.')], "_", str(i)]))

        if pep_df is None:
            pep_df = pd.DataFrame(data=zip(pep_data['pep_names'], pep_data['begInd'], pep_data['endInd'], peptides), columns=columns)
        else:
            pep_df = pd.concat((pep_df, pd.DataFrame(data=zip(name, pep_data['begInd'], pep_data['endInd'], peptides), columns=columns)), axis=0)


    # % indexStrains = [1:length(Seqs)];  
    # % for j=1:length(indexStrains)
    # %     currInds = strmatch(Seqs(indexStrains(j)).Header,{peptideStruct.Header});
        
    # %     if(j==1) 
    # %       inds = currInds;
    # %     else
    # %       currPeptides = {peptideStruct(currInds).sequence};
    # %       [C I] = setdiff(currPeptides,peptides);
    # %       inds = [inds; currInds(sort(I))]; 
    # %     end
    # %     peptides = {peptideStruct(inds).sequence};
        
    # % end
    # % peptideStruct = peptideStruct(inds);  

    if sort_method == 'beg_ind':
        pep_df.sort(columns='begInd', inplace=True)
    elif sort_method == 'header':   
        pep_df.sort(columns=['Header', 'begInd'], inplace=True)
    
    # output to file
    pep_df.to_csv(path_or_buf=output_filename, sep='\t', index=False)
    ipdb.set_trace()


    return pep_df
  




