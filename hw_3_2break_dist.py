from Bio import SeqIO
import argparse
from tqdm import tqdm

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in seq[::-1]])


def make_dotplot(seq1, seq2, k_size):
    
    # create dict with kmers indexes from seq1
    kmers_seq1_dict = {}
    # kmers_seq1_dict_rc = {}
    for s1_ind in tqdm(range(len(seq1) - k_size+ 1)):
        kmer = seq1[s1_ind:s1_ind+ k_size]

        if kmer not in kmers_seq1_dict:
            kmers_seq1_dict[kmer] = set()
        kmers_seq1_dict[kmer].add(s1_ind)

        # rc_kmer = reverse_complement(kmer)
        # if rc_kmer not in kmers_seq1_dict_rc:
        #     kmers_seq1_dict_rc[rc_kmer] = set()
        # kmers_seq1_dict_rc[rc_kmer].add(s1_ind)
    print('seq1 kmers created')
    
    dotplot_exact = []
    dotplot_rc = []

    for s2_ind in tqdm(range(len(seq2) - k_size+ 1)):
        kmer = seq2[s2_ind:s2_ind+k_size]

        # exact match: s2 kmer in s1 kmers
        if kmer in kmers_seq1_dict:

            for s1_ind in kmers_seq1_dict[kmer]:
                # dotplot_rc[kmer].add((s1_ind, s2_ind))
                dotplot_exact.append((s1_ind, s2_ind))
        
        # reverce complement (rc) match: s2 kmer in rc_s1 kmers or rc_s2 kmer in s1_kmers 
        reverse_complement_kmer = reverse_complement(kmer)
        if reverse_complement_kmer in kmers_seq1_dict:

            for s1_ind in kmers_seq1_dict[reverse_complement_kmer]:
                # dotplot_rc[kmer].add((s1_ind, s2_ind))
                dotplot_rc.append((s1_ind, s2_ind))

    return dotplot_exact, dotplot_rc


def remove_used_nodes(syntency_block, curr_nodes):
    for node in syntency_block:
        curr_nodes.remove(node)
    return 


def find_near_nodes(syntency_block, curr_nodes, max_distance):
    
    syntency_block.append(curr_nodes[0])
    
    for node in curr_nodes[1:]:
        print(curr_nodes)

        if node[0] - curr_nodes[0][0] > max_distance or node[1] - curr_nodes[0][1] > max_distance:
            continue
        if node[0] - curr_nodes[0][0] > max_distance and node[1] - curr_nodes[0][1] > max_distance:
            break
        print('add')
        syntency_block.append(node)
        find_near_nodes(syntency_block, curr_nodes[1:], max_distance)



def find_syntency_blocks(dotplot, max_distance, min_syntency_block_size):
        syntency_blocks_seq1 = []
        syntency_blocks_seq2 = []
        syntency_blocks = []
        syntency_block = [dotplot[0]]
        
        curr_node_ind= 0
        next_node_ind = 1
        # for next_node_ind in range(1, len(dotplot)):
        while next_node_ind < len(dotplot):
                
                # print(next_node_ind)
                node = dotplot[curr_node_ind]
                next_node = dotplot[next_node_ind]

                # print(node, next_node)
                # print('diffs',next_node[0]- node[0], next_node[1]- node[1])
                
                if abs(next_node[0] - node[0]) <= max_distance and abs(next_node[1] - node[1]) <= max_distance:
                        # print('add')
                        syntency_block.append(next_node)
                        curr_node_ind = next_node_ind

                if abs(next_node[0] - node[0]) > max_distance and abs(next_node[1] - node[1]) > max_distance:
                        
                        if len(syntency_block) >= min_syntency_block_size:
                                syntency_blocks.append(syntency_block)
                                syntency_blocks_seq1.append(min(syntency_block)[0])
                                syntency_blocks_seq2.append(min([(i[1], i[0]) for i in syntency_block])[0])
                                
                        curr_node_ind +=1
                        next_node_ind = curr_node_ind +1
                        syntency_block = [dotplot[curr_node_ind]]
                        continue

                next_node_ind +=1


        if len(syntency_block) >= min_syntency_block_size:
                syntency_blocks.append(syntency_block)
                syntency_blocks_seq1.append(min(syntency_block)[0])
                syntency_blocks_seq2.append(min([(i[1], i[0]) for i in syntency_block])[0])

        # return syntency_blocks
        return syntency_blocks, syntency_blocks_seq1, syntency_blocks_seq2


def make_graphs(
        syntency_blocks_seq1,
        syntency_blocks_rc_seq1,
        syntency_blocks_seq2,
        syntency_blocks_rc_seq2
):

    graph_seq1 = syntency_blocks_seq1 + syntency_blocks_rc_seq1
    graph_seq1 = [(value, index) for index, value in enumerate(graph_seq1)]
    graph_seq1_sorted = sorted(graph_seq1)
    graph_seq1_sorted = [node[1] for node in graph_seq1_sorted]


    graph_seq2 = syntency_blocks_seq2 + syntency_blocks_rc_seq2
    graph_seq2 = [(value, index) for index, value in enumerate(graph_seq2)]
    graph_seq2_sorted = sorted(graph_seq2)
    graph_seq2_sorted = [node[1] for node in graph_seq2_sorted]

    negative = [i for i in range(len(syntency_blocks_seq2), len(syntency_blocks_seq2) + len(syntency_blocks_rc_seq2))]
    return graph_seq1_sorted, graph_seq2_sorted, negative



def create_adj_list(graph_seq1_sorted, graph_seq2_sorted, negative):
    adj_dict = {}
    for node_ind, node in enumerate(graph_seq1_sorted):
        
        tmp_nodes = graph_seq1_sorted[-1:] + graph_seq1_sorted + graph_seq1_sorted[0:1]
        tmp_nodes = tmp_nodes[node_ind: node_ind+3]


        adj_dict[str(node) + 'in'] = [str(tmp_nodes[0]) + 'out']
        adj_dict[str(node) + 'out'] = [str(tmp_nodes[2]) + 'in']

    for node_ind, node in enumerate(graph_seq2_sorted):

        tmp_nodes = graph_seq2_sorted[-1:] + graph_seq2_sorted + graph_seq2_sorted[0:1]
        tmp_nodes = tmp_nodes[node_ind: node_ind+3]
        
        if node in negative:
            if tmp_nodes[0] in negative:
                adj_dict[str(node) + 'out'] += [str(tmp_nodes[0]) + 'in']
            else:
                adj_dict[str(node) + 'out'] += [str(tmp_nodes[0]) + 'out']
            
            if tmp_nodes[2] in negative:
                adj_dict[str(node) + 'in'] += [str(tmp_nodes[2]) + 'out']
            else:
                adj_dict[str(node) + 'in'] += [str(tmp_nodes[2]) + 'in']
        
        else:
            if tmp_nodes[0] in negative:
                adj_dict[str(node) + 'in'] += [str(tmp_nodes[0]) + 'in']
            else:
                adj_dict[str(node) + 'in'] += [str(tmp_nodes[0]) + 'out']
            
            if tmp_nodes[2] in negative:
                adj_dict[str(node) + 'out'] += [str(tmp_nodes[2]) + 'out']
            else:
                adj_dict[str(node) + 'out'] += [str(tmp_nodes[2]) + 'in']
    
    return adj_dict


def dfs(node, visited, adjacency_list):
    visited.add(node)
    for neighbor in adjacency_list[node]:
        if neighbor not in visited:
            dfs(neighbor, visited, adjacency_list)

def count_cycles(adjacency_list):
    visited = set()
    cycles = 0

    for node in adjacency_list:
        if node not in visited:
            dfs(node, visited, adjacency_list)
            cycles += 1

    return cycles

def load_data():
    human_seqs_path = 'data/chrX_human.fa'
    human_seqs = []
    for record in SeqIO.parse(human_seqs_path, "fasta"):
        human_seqs.append(str(record.seq))

    human_seq = human_seqs[0]
    human_seq = ''.join([i.upper() for i in human_seq if i.upper() != 'N'])

    mouse_seqs_path = 'data/chrX_mouse.fa'
    mouse_seqs = []
    for record in SeqIO.parse(mouse_seqs_path, "fasta"):
        mouse_seqs.append(str(record.seq))

    mouse_seq = mouse_seqs[0]
    mouse_seq = ''.join([i.upper() for i in mouse_seq if i.upper() != 'N'])
    return human_seq, mouse_seq


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-k', '--kmer_size', default=100)
    parser.add_argument('-d', '--max_distance', default=1)
    parser.add_argument('-b', '--min_syntency_block_size', default=5)
    parser.add_argument('-t', '--test', action="store_true")
    args = parser.parse_args()

    kmer_size = int(args.kmer_size)
    max_distance = int(args.max_distance)
    min_syntency_block_size = int(args.min_syntency_block_size)

    print('0: start loading data')
    if args.test:
        mouse_seq = 'AGCAGGAGATAAACCTGT'
        human_seq = 'AGCAGGTTATCTACCTGT'
    else:
        human_seq, mouse_seq = load_data()
    print('\tlen seq1:', len(human_seq))
    print('\tlen seq2:', len(mouse_seq))

    print('1: start making dotplot')
    dotplot, dotplot_rc = make_dotplot(mouse_seq, human_seq, kmer_size)
    dotplot_sorted = sorted(dotplot)
    dotplot_rc_sorted = sorted(dotplot_rc)
    print('\tnumber of points in dotplot forward:', len(dotplot))
    print('\tnumber of points in dotplot reverse:', len(dotplot_rc))

    print('2: start making syntency blocks')
    syntency_blocks, syntency_blocks_seq1, syntency_blocks_seq2 = find_syntency_blocks(dotplot_sorted, max_distance, min_syntency_block_size)
    syntency_blocks_rc, syntency_blocks_rc_seq1, syntency_blocks_rc_seq2 = find_syntency_blocks(dotplot_rc_sorted, max_distance, min_syntency_block_size)
    print('\tnumber of syntency_blocks forward:', len(syntency_blocks))
    print('\tnumber of syntency_blocks reverse:', len(syntency_blocks_rc))

    print('3: start making graphs')
    graph_seq1_sorted, graph_seq2_sorted, negative = make_graphs(
            syntency_blocks_seq1,
            syntency_blocks_rc_seq1,
            syntency_blocks_seq2,
            syntency_blocks_rc_seq2)

    print('4: start making adj list')
    adj_dict = create_adj_list(graph_seq1_sorted, graph_seq2_sorted, negative)
    print('5: start counting cycles')
    cycles = count_cycles(adj_dict)

    print('***number of syntency blocks', len(graph_seq1_sorted))
    print('***number of cycles', cycles)
    print('***FINAL RESULT', len(graph_seq1_sorted) - cycles)


