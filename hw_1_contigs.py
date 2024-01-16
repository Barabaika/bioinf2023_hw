from tqdm import tqdm
import argparse


'''
Script for predicting contigs from two fasta files.
Need tqdm package.

ARGS:
1- first input fasta
2- second input fasta

KWARGS:
-k --kmer_size (default = 51) : size of created kmers from reads

NOTES:
- very slow (took more than 24 hours (about 1500 minutes) to process hw1 files)

USAGE:
python hw_1_contigs.py data/Carsonella_ruddii_reads_paired_reads_left.fastq data/Carsonella_ruddii_reads_paired_reads_right.fastq -o hw1_test_out.txt -k 51

'''

def get_kmer_count_from_sequence(sequence, k=3, cyclic=True):
    """
    Returns dictionary with keys representing all possible kmers in a sequence
    and values counting their occurrence in the sequence.
    """
    # dict to store kmers
    kmers = {}
    
    # count how many times each occurred in this sequence (treated as cyclic)
    for i in range(0, len(sequence)):
        kmer = sequence[i:i + k]
        
        # for cyclic sequence get kmers that wrap from end to beginning
        length = len(kmer)
        if cyclic:
            if len(kmer) != k:
                kmer += sequence[:(k - length)]
        
        # if not cyclic then skip kmers at end of sequence
        else:
            if len(kmer) != k:
                continue
        
        # count occurrence of this kmer in sequence
        if kmer in kmers:
            kmers[kmer] += 1
        else:
            kmers[kmer] = 1
    
    return kmers


def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in seq[::-1]])

def add_tmp_dict_to_res_dict_kmer(tmp_kmers, all_kmers):
    for key in tmp_kmers.keys():
        # if there is already this kmer in all dict - add count to its value
        if key in all_kmers.keys():
            all_kmers[key] += tmp_kmers[key]
        else:
            all_kmers[key] = tmp_kmers[key]
    return all_kmers


def get_kmer_count_from_paired_fasta(reads_file1, reads_file2, k):
    
    all_kmers = {}
    with open(reads_file1) as ff, open(reads_file2) as fr:
        for linef, liner in zip(ff, fr):
            if linef.startswith('@'):  # Assuming standard Fastq format
                sequence1 = next(ff).strip()
                sequence2 = next(fr).strip()
                
                # # reverce complement to be able to use both forward and reverse reads for assembling one forward strand
                # sequence_reverse_f = reverse_complement(sequence_reverse)
                kmers1 = get_kmer_count_from_sequence(sequence1, k, cyclic=False)
                kmers2 = get_kmer_count_from_sequence(sequence2, k, cyclic=False)

                # store all kmers in result dict (all_kmers) 
                all_kmers = add_tmp_dict_to_res_dict_kmer(kmers1, all_kmers)
                all_kmers = add_tmp_dict_to_res_dict_kmer(kmers2, all_kmers)                

    return all_kmers


def get_debruijn_nodes_from_kmers(kmers):

    # store nodes in dict, in format dict[paremt_node] = child_node
    nodes = {}
    
    for k in kmers:
        if k[:-1] not in nodes.keys():
            nodes[k[:-1]] = set()
        nodes[k[:-1]].add(k[1:])
    return nodes


def count_in_nodes(node_name, nodes):
    counter= 0
    for node in nodes:
        if node_name in nodes[node]:
           counter +=1 
    return counter

def is_it_one_in_one_out_node(node_name, nodes):
    in_number = count_in_nodes(node_name, nodes)
    if node_name in nodes:
        out_number = len(nodes[node_name]) 
    else:
        out_number = 0
    return (in_number == 1 and out_number == 1)


def get_one_element_from_set(node_set):
    for e in node_set:
        break
    return e


def maximal_nonbranching_paths(nodes):
    paths = []
    for node in tqdm(nodes):
        # print(node)
        # print(nodes[node])
        out_number = len(nodes[node])
        # check if element is one in - one out:
        if not is_it_one_in_one_out_node(node, nodes):
            if out_number > 0:
                for child_node in nodes[node]:
                    nonbranchingpath = node + child_node[-1]
                    while is_it_one_in_one_out_node(child_node, nodes):
                        
                        child_node = get_one_element_from_set(nodes[child_node])
                        nonbranchingpath += child_node[-1]
                    paths.append(nonbranchingpath)

    return paths


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("input_file1")
    parser.add_argument("input_file2")
    parser.add_argument('-o', '--outfile', default= 'out_hw1.txt')
    parser.add_argument('-k', '--kmer_size', default=51)
    args = parser.parse_args()

    print(f'start process {args.input_file1} and {args.input_file2}\n')

    path1 = args.input_file1 #'data/Carsonella_ruddii_reads_paired_reads_left.fastq'
    path2 = args.input_file2 #'data/Carsonella_ruddii_reads_paired_reads_right.fastq'
    
    # create kmers from reads in fasta files 
    all_kmers = get_kmer_count_from_paired_fasta(path1, path2, k = int(args.kmer_size))
    # create dict with nodes from kmers, this is not a full debruijn graph. 
    # May be need to create exactly debruijn graph to reduce time usage
    all_nodes = get_debruijn_nodes_from_kmers(all_kmers)
    # get maximal nonbranching paths in debruijn graph == contigs
    paths = maximal_nonbranching_paths(all_nodes)
    # store contigs to txt file
    with open(args.outfile, 'w') as f:
        for p in paths:
            f.write(p + '\n')

    print(f'Done, number of contigs: {len(paths)}')


