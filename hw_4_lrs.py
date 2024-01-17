import argparse
'''
Script to find longest repeted substring (LRS)

ARGS:
1- input sequense: string

NOTES:
- LRS can be overlaped
- no output file, only print in consol  

USAGE:
python hw_4_lrs.py banana
'''

def create_suffix_array(seq):
    suffix_array = [(seq[i:], i) for i in range(len(seq))]
    suffix_array.sort()
    suffix_array = [i[1] for i in suffix_array]
    return suffix_array


def create_lcp(suffix_array, seq):
    # rank is a positin of suffix in suffix array. ex: rank[0] is a position of the first suffix (full sequence) in sorted suffix array
    # with rank we can start from the begining of the sequence and know what is the index if this suffix      
    rank = [0]* len(seq)
    for i in range(len(seq)):
        rank[suffix_array[i]] = i
    
    lcp_values = [0]*len(seq)
    h= 0 # is a lcp length
    # go from full sequence (suffix_array[rank[0]]) to the previous to last suffix
    for i in range(len(seq)- 1):
        prev_suffix_ind = suffix_array[rank[i]-1]
        while i+h < len(seq) and prev_suffix_ind + h < len(seq) and seq[prev_suffix_ind+h] == seq[i+h]:
            h+=1
        
        lcp_values[rank[i]] = h
        # we know that the next suffix can have min lcp-1, so for h next we need only to do h- 1
        if h >0:
            h -=1
    
    return lcp_values


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("input_seq")
    args = parser.parse_args()

    print(f'start process\nsequence length = {len(args.input_seq)}')

    seq= args.input_seq + '$'
    suffix_array = create_suffix_array(seq)
    lcp_array = create_lcp(suffix_array, seq)

    max_length_lcp = max(lcp_array)
    max_index = lcp_array.index(max_length_lcp)

    lrs = seq[suffix_array[max_index]:suffix_array[max_index]+max_length_lcp]
    print(f'Longest Repeated Substring: {lrs}')



