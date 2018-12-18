'''
synbiochem (c) University of Manchester 2018

synbiochem is licensed under the MIT License.

To view a copy of this license, visit <http://opensource.org/licenses/MIT/>.

@author:  neilswainston
'''
from collections import defaultdict
import json
from random import choice
import sys

from Bio import SeqIO
from Bio.Alphabet import generic_protein
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import matplotlib.pyplot as plt
import pandas as pd
from synbiochem.utils import net_utils


def get_fasta(prefix, suffix, entries, out_file='out.fasta',
              organism='4932', seqs_per_class=100):
    '''Get fasta file of protein variants.'''
    base_url = 'http://codon.synbiochem.co.uk/codons?codon=%s&organism=%s'

    class_seqs = defaultdict(list)

    for entry in entries:
        seq_class, nucl_seq = entry.split('_')
        whole_seq = prefix + nucl_seq + suffix
        num_of_codons = len(whole_seq) // 3
        codon_indexes = range(0, num_of_codons * 3, 3)

        aa_seqs = []

        for codon_index in codon_indexes:
            codon = whole_seq[codon_index:codon_index + 3]
            url = base_url % (codon, organism)
            resp = json.loads(net_utils.get(url, headers=None,
                                            verify=True))
            aas = [aa['amino_acid'] for aa in resp[0]['amino_acids']
                   if aa['amino_acid'] != 'Stop']

            if aa_seqs:
                new_aa_seq = []

                for aa in aas:
                    for aa_seq in aa_seqs:
                        new_aa_seq.append(aa_seq + aa)

                aa_seqs = new_aa_seq
            else:
                aa_seqs = aas

        class_seqs[seq_class].extend(aa_seqs)

    with open(out_file, 'w') as out:
        for seq_class, seqs in class_seqs.items():
            seqs = [choice(seqs)
                    for _ in range(min(len(seqs), seqs_per_class))]

            SeqIO.write([SeqRecord(Seq(aa_seq, generic_protein),
                                   '%s_%i' % (seq_class, aa_idx + 1),
                                   '', '')
                         for aa_idx, aa_seq in enumerate(seqs)],
                        out, 'fasta')


def analyse(filename, valid_range):
    '''Analyse results.'''
    columns = ['Class assignment - B for buried or E for Exposed - ' +
               'Threshold: 25% exposure, but not based on RSA',
               'Amino acid',
               'Sequence name',
               'Amino acid number',
               'Relative Surface Accessibility - RSA',
               'Absolute Surface Accessibility',
               'Z-fit score for RSA prediction',
               'Probability for Alpha-Helix',
               'Probability for Beta-strand',
               'Probability for Coil']

    df = pd.read_csv(filename, comment='#', sep='\s+')
    df.columns = columns

    names = df['Sequence name'].str.split('_', n=2, expand=True)
    df['class'] = names[0]
    df['variant_index'] = names[1]

    valid_aas = list(sorted(df['Amino acid number'].unique()))[
        valid_range.start:valid_range.stop]

    valid_df = df[df['Amino acid number'].isin(valid_aas)]

    class_aa_nums = []

    for idxs, group_df in valid_df.groupby(['class', 'Amino acid number']):
        class_aa_num = group_df['Probability for Alpha-Helix'].describe()
        class_aa_num['class'] = idxs[0]
        class_aa_num['Amino acid number'] = idxs[1] - valid_range.start
        class_aa_nums.append(class_aa_num)

    summary_df = pd.concat(class_aa_nums, axis=1).T
    summary_df.to_csv('summary.csv')
    return summary_df


def plot(df, out_filename='struct_preds.png'):
    '''Plot.'''
    plt.rcParams['font.family'] = 'Times New Roman'

    for clss, group_df in df.groupby('class'):
        plt.errorbar(group_df['Amino acid number'], group_df['mean'],
                     yerr=group_df['std'],
                     label=clss,
                     fmt='o', alpha=0.5, capsize=3)

    plt.title('Probability for Alpha-Helix predictions')
    plt.xlabel('Amino acid number')
    plt.ylabel('Probability for Alpha-Helix')
    plt.legend(loc='upper right')
    plt.savefig(out_filename, bbox_inches='tight')


def main(args):
    '''main method.'''
    # get_fasta(args[0], args[1], args[2:])
    summary_df = analyse('data/result.txt',
                         range(len(args[0]) // 3, len(args[1]) // -3))
    plot(summary_df)


if __name__ == '__main__':
    main(sys.argv[1:])
