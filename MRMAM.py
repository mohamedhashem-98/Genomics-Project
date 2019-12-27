def readFastq(filename):
    sequences = []
    qualities = []
    fh = open(filename)
    seq_length = fh.readline()
    while True:
        fh.readline()  # skip ID line
        seq = fh.readline().rstrip()  # read base sequence
        fh.readline()  # skip + sign
        fh.readline().rstrip()  # base quality line
        if len(seq) == 0:
            break
        sequences.append(seq)
    return seq_length, sequences


def readtextfile(filename):
    length = 0
    f = open(filename, "r")
    sequences = []
    length = f.readline().split('\n')[0]
    while True:
        seq = f.readline().split('\n')[0]
        if len(seq) == 0:
            break
        sequences.append(seq)
    return length, sequences


def new_kmers(kmers):   # This function generates the K-1-Mers, given the K-mers.
    length = len(kmers[0])-1
    new_list = []
    val = 0
    tmp = []
    my_dict = {}
    count = 1
    for kmer in kmers:
        k = count
        tmp = []
        val = 0
        while True:
            new_kmer = kmer[val:val+length]     # Generating K-1-Mers to every Kmer.
            if len(new_kmer) != length:
                break
            else:
                tmp.append(new_kmer)
            val += 1
        new_list.append(tmp)
        my_dict[k] = tmp            # Storing each 2 K-1-mers in a dictionary with 1-based index as key.
        count += 1
    return my_dict


# This function builds the Graph given the dictionary of K-1-mers.
def graph(my_dict):
    nodes = []
    nodes.append(my_dict[1][0])
    edges = []
    new_dic = {}
    for key, value in my_dict.items():      # Each node is distinct.
        for val in value:
            if val in nodes:
                continue
            else:
                nodes.append(val)

    for key, value in my_dict.items():      # Each value in the given dictionary represents the edges of the graph.
        edges.append(value)

    tmp = []
    for nod in nodes:
        tmp = []
        for li in edges:
            if nod == li[0]:        # Here, we append the neighbors of each nod
                tmp.append(li[1])
        if len(tmp) == 0:
            tmp.append(-1)          # Here, we determine the end node of the graph by putting -1 in the value of this node.
        new_dic[nod] = tmp          # We store each node as a key, with it's neighbors as value in form of list.

    return new_dic


def Eulerian(graph):
    my_keys = list(graph.keys())
    my_values = list(graph.values())
    new_list = []
    path = ''

    for li in my_values:
        for i in range(len(li)):
            new_list.append(li[i])

    for k in my_keys:           # Here, we determine the start of the graph, by finding which node is not found in any
        if k not in new_list:   # adjacent list of any other node.
            start = k

    path += start
    li = []
    for key, value in graph.items():    # By traversing the graph, we pick the next node which is adjacent to the start node                         #$
        for i in range(len(value)):     # and the start node is updated along the loop.
            start = graph[start][i]     # Here, we're handling the case if the node has adjacency list of length > 1
            if start == -1:
                break                   # if we reached the end node (-1) the loop breaks.
            li.append(start)            # Appending each node in order in a temp list.

    for i in range(len(li)):
        path += li[i][-1]               # Taking the last character of each node.

    return path


# Multiple reads:

def readtextfilepaired(filename):
    lengthofseq = 0
    f = open(filename, "r")
    sequances = []
    firstline = f.readline().split(' ')
    lengthofseq = int(firstline[0])
    lengthofgap = int(firstline[1])
    suffixlength = lengthofseq + lengthofgap
    while True:
        seq = f.readline().split('\n')[0]
        if len(seq) == 0:
            break
        sequances.append(seq)
    return sequances, suffixlength


def new_kmerspaired(sequance):
    tmp = ''
    before = ''
    after = ''
    dic = {}
    for i in range(len(sequance)):
        tmp = sequance[i]
        before = tmp.split('|')[0]    #before contain sequence before |
        after = tmp.split('|')[1]     #after contain sequence after |
        tmp1 = []
        tmp2 = []
        for val in range(len(before) - 1):
            new_before = before[val:val + len(before) - 1]      #generate k-1Mers of sequence before |
            new_after = after[val:val + len(after) - 1]         #generate k-1Mers of sequence after |
            tmp1.append(new_before)                             #appending k-1Mers of before in tmp1
            tmp2.append(new_after)                              #appending k-1Mers of after in tmp2

        dic[tmp1[0], tmp2[0]] = (tmp1[1], tmp2[1])              #key of dic is the prefix, value is the suffix
    return dic


def Eulirean_Paired(mydic, s):
    startp = ()
    starts = ()
    suffix = []
    prefix = []
    prefixpath = []
    suff_path = []
    prefix_str = ''
    suffix_str = ''
    path = ''
    for key, value in mydic.items():
        prefix.append(key)
        suffix.append(value)
    for i in range(len(prefix)):                #determine the start prefix and start suffix that is not found in adjacent list of any node
        if prefix[i] not in suffix:
            startp = prefix[i]
            starts = prefix[i]
    prefixpath.append(startp[0])
    suff_path.append(starts[1])
    for key, value in mydic.items():            #by traversing the grapth we update the start node that is adjacent of the current start
        startp = mydic[startp]
        starts = mydic[starts]
        prefixpath.append(startp[0])
        suff_path.append(starts[1])
    for i in range(len(prefixpath) - 1):
        prefix_str += prefixpath[i][0]                      # Appending the prefix nodes in order
        suffix_str += suff_path[i][0]                       # Appending the suffix nodes in order
    prefix_str += prefixpath[len(prefixpath) - 1]           # Taking last character of each prefix
    suffix_str += suff_path[len(suff_path) - 1]
    path = prefix_str + suffix_str[len(suffix_str) - s::]   # After generating the suffix, prefix we put the prefix with the last k+d of suffix
    return path


while True:
    ch = input('Single or Multiple read         [1: Single   2: Multiple   0: End] \n')
    if ch == '1':
        length_of_seq, sequence = readtextfile("input.txt")
        new_km = new_kmers(sequence)
        G = graph(new_km)
        p = Eulerian(G)
        f = open('MyOutput.txt', 'w')
        f.write(p)
        f.close()
        print('Output is exported to .txt file \n')
    elif ch == '2':
        sequance, suffixlength = readtextfilepaired("Paired_Input.txt")
        # sequance=['GACC|GCGC','ACCG|CGCC','CCGA|GCCG','CGAG|CCGG','GAGC|CGGA']
        D = new_kmerspaired(sequance)
        p = Eulirean_Paired(D, suffixlength)
        f = open("PairedOutput.txt", "w")
        f.write(p)
        f.close()
        print('Output is exported to .txt file \n')
    else:
        print('Thank you \n')
        break








