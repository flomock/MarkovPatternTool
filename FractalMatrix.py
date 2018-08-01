import re
import math
import matplotlib

# Force matplotlib to not use any Xwindows backend.
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from os import listdir
from os.path import isfile, join
import editdistance
from multiprocessing import Pool
import multiprocessing
from itertools import repeat

c_score = 0
g_score = 0
t_score = 0
gc_score = 0
plt.rcParams['figure.figsize'] = 10, 10


def one_hot_encoding(seq):
    """
    1. parses seq to binary coordinates
    2. calculate the C,G,A,T content
    :param seq: input genome sequence
    :return: binary sequence
    """

    global g_score
    global c_score
    global t_score

    seq_array = np.array(list(seq))
    x_bin = np.zeros(seq_array.shape, dtype=np.int8)
    y_bin = np.zeros(seq_array.shape, dtype=np.int8)
    z_bin = np.zeros(seq_array.shape, dtype=np.int8)

    y_bin[seq_array == 'A'] = 1
    a_score = y_bin.sum()
    y_bin[seq_array == 'G'] = 1
    g_score = y_bin.sum() - a_score

    z_bin[seq_array == 'A'] = 1
    z_bin[seq_array == 'C'] = 1
    c_score = z_bin.sum() - a_score
    t_score = len(seq) - a_score - c_score - g_score

    return x_bin, y_bin, z_bin


def buildMatrix(size):
    """
    build matrix, each cell representing one tuple / k-mer
    :param size: k-mer length
    :return: empty matrix
    """
    # global matrix
    matrix = np.zeros((2 ** size, 2 ** size), dtype=np.int)
    return matrix


def readDNA(path, endsize, subseq=False, filter_mikroSats=False):
    """
    1. read the input file
    2. prepare big matrix and parsed sequence
    3. parse sequence to entry in matrix
    :param path: path to input file
    :param endsize: size of tuple
    :param subseq: set True if only fraction of sequence is needed
    :return: matrix with entries of tuple
    """
    global gc_score
    global c_score
    global g_score
    global t_score
    gc_score = 0

    matrix = buildMatrix(endsize)

    if subseq == False:

        with open(path) as f:
            dnaSeq = ""
            iterlines = iter(f)
            first_line = next(iterlines)
            if first_line.startswith(">"):
                pass
            else:
                dnaSeq += first_line.upper()[:-1]

            for line in f:
                dnaSeq += line.upper()[:-1]

            subseqs_DNA = re.split(r"[^ACGT]+", dnaSeq)

    else:
        subseqs_DNA = re.split(r"[^ACGT]+", subseq)

    numCPUs = multiprocessing.cpu_count()
    pool = Pool(numCPUs)
    if filter_mikroSats:
        subseqs_DNA = pool.map(del_microsats, subseqs_DNA)
        subseqs_DNA = [x for y in subseqs_DNA for x in y]

    # print(''.join(subseqs_DNA))
    # exit()

    x_bins, y_bins, z_bins = [], [], []
    for seq in subseqs_DNA:
        x_binary, y_binary, z_binary = one_hot_encoding(seq)
        x_bins.append(x_binary)
        y_bins.append(y_binary)
        z_bins.append(z_binary)

    input_bins = []
    for cpu in range(numCPUs):
        input_bins.append([])

    for index, i in enumerate(range(len(y_bins))):
        input_bins[index % numCPUs].append([y_bins[i], z_bins[i]])

    matrix_list = pool.starmap(countTuple, zip(input_bins, repeat(endsize)))

    matrix = np.sum(matrix_list, axis=0)
    totalLength = np.sum(matrix)

    c_score = c_score / totalLength
    g_score = g_score / totalLength
    t_score = t_score / totalLength

    return matrix


def bin_to_int(binary, biggest_bin, old=-1):
    """
    parse binary numbers to int and use fact that old already calculated
    :param binary: input sequence of binary numbers
    :param biggest_bin: biggest possible number
    :param old: old decimal integer number
    :return: decimal integer
    """
    decimal = old % biggest_bin
    decimal = decimal * 2 + int(binary[-1])
    if (old == -1):
        decimal = 0
        for digit in binary:
            decimal = decimal * 2 + int(digit)
    return decimal


def countTuple(input, endsize):
    matrix_sub = buildMatrix(endsize)
    for i in range(len(input)):
        biggest_bin = 2 ** (endsize - 1)
        old_y = -1
        old_z = -1
        input_bin = input[i]
        y_binary = input_bin[0]
        z_binary = input_bin[1]

        for posi in range(0, len(y_binary) - endsize + 1, 1):
            y = bin_to_int(y_binary[posi:posi + endsize], biggest_bin, old_y)
            z = bin_to_int(z_binary[posi:posi + endsize], biggest_bin, old_z)
            old_y = y
            old_z = z

            matrix_sub[y][z] += 1

    return matrix_sub


# ACHTUNG gilt nur fuer num i = y und num j = z
def int_to_binarray(num_i, num_j, size):
    """
    parse decimal to binary number
    :param num_i: first number to parse (y axis)
    :param num_j: second number to parse (z axis)
    :param size: tuple size
    :return: binary numbers
    """
    binary_i = "{0:b}".format(num_i)
    binary_j = "{0:b}".format(num_j)
    foo = np.zeros((2, size), dtype=np.int)
    # array mit fuehrenden nullen erstellen
    binary_i = binary_i[::-1]
    binary_j = binary_j[::-1]

    for i in range(0, len(binary_i), 1):
        foo[0, i] = binary_i[i]
    for i in range(0, len(binary_j), 1):
        foo[1, i] = binary_j[i]

    foo[0] = foo[0, ::-1]
    foo[1] = foo[1, ::-1]
    return foo


def plot_image(outMatrix, outputName, sequence_size=None):
    """
    plot the matrix as pdf
    :param outMatrix: matrix to be plotted
    :param outputName: name of the to be created file
    :param sequence_size: length of the sequence
    :return: -
    """
    if sequence_size == None:
        sequence_size = int((round(np.sum(outMatrix) / 100)) * 100)
    m = outMatrix.flatten()
    m_sorted = np.sort(m)
    nonZeroSorted = m_sorted[m_sorted.nonzero()]
    opaque = nonZeroSorted[int(0.95 * len(nonZeroSorted))]  # cut off bei 95% von nicht leeren Feldern
    pic = [[float(outMatrix[row][col]) / float(opaque) for row in range(len(outMatrix[0]))] for col in
           range(len(outMatrix[0]))]
    for i in range(len(pic)):
        for j in range(len(pic)):
            if pic[i][j] > 1:
                pic[i][j] = 1

    plt.imshow(pic, cmap='gray_r', interpolation="none", extent=[0, 1, 1, 0])
    # show which character belongs to which direction
    # plt.text(-0.08,1.03,"C",fontsize=35)
    # plt.text(1.03,1.03,"A",fontsize=35)
    # plt.text(-0.08,-0.08,"T",fontsize=35)
    # plt.text(1.03,-0.08,"G",fontsize=35)
    plt.gca().invert_yaxis()
    plt.title("C:" + str(int(round(c_score * 100))) + "%  " + "G:" + str(int(round(g_score * 100))) + "%  " +
              "T:" + str(int(round(t_score * 100))) + "%  " + "A:" + str(
        int(round((1 - t_score - c_score - g_score) * 100))) + "%\n" +
              "ca. " + str(sequence_size) + "bp")

    ax = plt.gca()
    ax.set_autoscale_on(False)

    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(outputName + ".pdf")
    plt.close()


def int_to_seq(num_i, num_j, size):
    """
    parses an integer combination to the resembling sequence
    :param num_i: first number to parse (y axis)
    :param num_j: second number to parse (z axis)
    :param size: tuple size
    :return: sequence text
    """
    foo = int_to_binarray(num_i, num_j, size)
    seq = []
    for k in range(0, size):
        if foo[0, k] == 0 and foo[1, k] == 0:
            seq.append("T")
        elif foo[0, k] == 0 and foo[1, k] == 1:
            seq.append("C")
        elif foo[0, k] == 1 and foo[1, k] == 0:
            seq.append("G")
        elif foo[0, k] == 1 and foo[1, k] == 1:
            seq.append("A")
    return seq


def sorting(outMatrix, size, outputName, plotting, log):
    """
    sorting the tuple by the number of occurrences
    :param outMatrix: matrix to sort
    :param size: size of the tuples
    :param outputName: output-file-name
    :param plotting: True if plotting wanted
    :param log: True if log scale wanted
    :return: sorted list
    """

    sorted_seq = []
    m = outMatrix.flatten()
    m_sorted = np.argsort(m)
    if plotting == True:
        m_small = m_sorted[0:15]
        m_small = np.append(m_small, m_sorted[len(m_sorted) - 15:len(m_sorted)])
        m_sorted = m_small

    # calculate placement in array
    for number in m_sorted:
        column = number % (2 ** size)
        row = int((number - column) / (2 ** size))
        seq = [m[number], int_to_seq(row, column, size)]
        sorted_seq.append(seq)

    if plotting == True:
        if m.max() == 0:
            return

        # building plot
        ax = plt.subplot()

        ind = np.arange(len(m_small))
        name = []
        quantity = []
        for i in range(len(sorted_seq)):
            name.append("".join(sorted_seq[i][1]))
            quantity.append(sorted_seq[i][0])

        # plotting
        plt.yscale('log')
        ax.bar(ind, quantity, 0.2)
        plt.xticks(ind, name, rotation='vertical')

        plt.tight_layout()
        plt.ylim(ymin=1)
        if "/" not in outputName:
            plt.savefig('/home/go96bix/Dropbox/hiwiManja/Fraktale/FractalDNA/MostLeast/' + outputName + str(
                size) + "MostLeast" + ".png")
        else:

            plt.savefig(outputName + "-MostLeast" + ".png")
        plt.close()
    elif (plotting == "energy"):
        return sorted_seq
    else:
        if log:
            for i in range(0, len(sorted_seq)):
                sorted_seq[i][0] = np.round(np.log(sorted_seq[i][0]), decimals=4)
        if "/" not in outputName:
            with open('/home/go96bix/Dropbox/hiwiManja/Fraktale/FractalDNA/MostLeast/' + outputName + ".txt",
                      "w") as text_file:
                for item in sorted_seq:
                    text_file.write("%s" % item[0] + "\t" + str("".join(item[1])) + "\n")
        else:
            with open(outputName + ".txt", "w") as text_file:
                for item in sorted_seq:
                    text_file.write("%s" % item[0] + "\t" + str("".join(item[1])) + "\n")


def matrix_to_csv(outMatrix, StandartOutputName):
    """
    save matrix in csv file
    :param outMatrix: matrix to be saved
    :param StandartOutputName: file name or True
    :return: -
    """
    r = re.compile(r"/")
    outputName = r.sub("-", data) + ".csv"
    if StandartOutputName != True:
        # outputName = r.sub("-", data) + str(StandartOutputName) + ".csv"
        outputName = str(StandartOutputName) + ".csv"
    np.savetxt(outputName, outMatrix, fmt='%.1d', delimiter=',')


def csv_to_matrix(file, scoring):
    """
    load csv file as integer numpy matrix
    :param file: file path
    :param scoring: if True estimate GCAT content
    :return: matrix
    """
    print(file)
    newMatrix = np.genfromtxt(file, delimiter=',', dtype='int32')
    global c_score, g_score, t_score, gc_score
    if scoring == True:
        seq_length = 0
        word_length = int(math.log(len(newMatrix[0]), 2))
        global gc_score
        for i in range(0, len(newMatrix[0])):
            for j in range(0, len(newMatrix[0])):
                # gc_in_word = 0
                c_in_word = 0
                g_in_word = 0
                t_in_word = 0
                seq = int_to_seq(i, j, word_length)
                for foo in seq:
                    if foo == "C" or foo == "c":
                        c_in_word += 1
                    if foo == "G" or foo == "g":
                        g_in_word += 1
                    if foo == "T" or foo == "t":
                        t_in_word += 1
                c_score += c_in_word * newMatrix[i][j]
                g_score += g_in_word * newMatrix[i][j]
                t_score += t_in_word * newMatrix[i][j]
                # gc_score += gc_in_word * newMatrix[i][j]
                seq_length += word_length * newMatrix[i][j]
        if seq_length != 0:
            g_score = float(float(g_score) / float(seq_length))
            c_score = float(float(c_score) / float(seq_length))
            t_score = float(float(t_score) / float(seq_length))
            gc_score = g_score + c_score
        else:
            t_score = 0
            c_score = 0
            g_score = 0

    return newMatrix


def csvs_to_plot(mypath):
    """
    plot all csv's in path
    :param mypath: path to dir of csvs
    :return:
    """
    paths = [mypath + '/' + f for f in listdir(mypath) if isfile(join(mypath, f))]
    for i in paths:
        foo = csv_to_matrix(i, True)
        name = str(i)[len(mypath) + 1:-4]
        plot_image(foo, name)


def csvs_to_moLeCommon(mypath, size, plotting):
    """
    get most extrem outliers of multiple csv's
    :param mypath: path to csv files
    :param size: tuple size
    :param plotting: True if plotting wanted
    :return:
    """
    for root, dirs, files in os.walk(mypath):
        for filename in [f for f in files if f.endswith(".csv")]:
            address = str(os.path.join(root, filename))
            foo = csv_to_matrix(address, False)
            if plotting:
                if not os.path.isdir(str(root + "/MostLeast/")):
                    os.makedirs(str(root + "/MostLeast/"))
                sorting(foo, size, str(root + "/MostLeast/" + filename[:-4]), plotting)
            else:
                sorting(foo, size, address[:-4], plotting)


def csv_representation_level_markov_chain(input_matrix, length_k=2):
    """
    main method of pattern prediction with markov chains
    1. determine the probabilities for tuple appending
    e.g. how likely ACT --> ACTG
    2. make prediction how matrix with length n will look like
    3. return real- / predicted- data

    :param input_matrix: matrix to be analyzed
    :param length_k: length of k parameter
    :return: expression level matrix
    """

    size = int(math.log(len(input_matrix[0]), 2))
    matrix_k = matrix_shrink_size(length_k, input_matrix)
    matrix_k_minus1 = matrix_shrink_size(length_k - 1, input_matrix)
    #
    # transitionmatrix = determine the probabilities for tuple appending
    transition_matrix = np.zeros((2 ** length_k, 2 ** length_k), dtype=np.float)
    for i in range(0, 2 ** length_k):
        for j in range(0, 2 ** length_k):
            k = float(matrix_k[i, j])
            k_minus = float(matrix_k_minus1[int(i / 2), int(j / 2)])
            if k_minus == 0:
                transition_matrix[i, j] = 1
            else:
                transition_matrix[i, j] = float(k / k_minus)

    # make prediction how matrix with length n will look like
    prediction_matrix = np.divide(np.array(matrix_k).astype(float), np.sum(matrix_k) - length_k + 1)
    for i in range(1, size - length_k + 1):
        helper_matrix = np.zeros((2 ** (length_k + i), 2 ** (length_k + i)), dtype=np.float)

        for l in range(0, 2 ** (length_k + i - 1)):
            for m in range(0, 2 ** (length_k + i - 1)):
                for x in (0, 1):
                    for y in (0, 1):
                        # cut first letter (search in matrix at the right location)
                        tr_x = (l % 2 ** (length_k - 1)) * 2
                        tr_y = (m % 2 ** (length_k - 1)) * 2
                        prediction = prediction_matrix[l, m] * transition_matrix[tr_x + x, tr_y + y]
                        helper_matrix[l * 2 + x][m * 2 + y] = prediction
        prediction_matrix = helper_matrix
    # parse number of observations of tuples into how likelihood to see certain pattern
    input_matrix_prob = np.divide(input_matrix.astype(float), np.sum(input_matrix))
    # return real- / predicted- data
    out_matrix = np.divide(input_matrix_prob, prediction_matrix)
    # Replace nan with zero and inf with finite numbers.
    out_matrix = np.nan_to_num(out_matrix)

    return out_matrix


def matrix_shrink_size(new_size, input_matrix):
    """
    shrnk the input matrix to the dedicated size
    :param new_size: size wanted
    :param input_matrix: matrix to parse
    :return: new smaller matrix
    """
    size = int(math.log(len(input_matrix[0]), 2))
    out_matrix = buildMatrix(new_size)

    for i in range(0, 2 ** new_size):
        for j in range(0, 2 ** new_size):
            out = 0
            for x in range(0, (2 ** (size - new_size))):
                for y in range(0, (2 ** (size - new_size))):
                    out += input_matrix[i * (2 ** (size - new_size)) + x][j * (2 ** (size - new_size)) + y]
            out_matrix[i, j] = out
    return out_matrix


def path_to_markovPatternAnalyse(mypath, length_n, length_k, recursiv, log):
    """
    analyze multiple csv for inhibition and activation
    :param mypath: path to to dir or file
    :param length_n: length of the tuple's to analyze
    :param length_k: length of the tuple's to extend
    :param recursiv: if True analyze also sub-dirs
    :param log: if True use log scale
    :return:
    """

    if (os.path.isdir(mypath)):
        for root, dirs, files in os.walk(mypath):
            for filename in [f for f in files if f.endswith(".csv")]:
                address = str(os.path.join(root, filename))
                classical_matrix = csv_to_matrix(address, True)
                markovPatternAnalyse(path=address, classical_matrix=classical_matrix, length_n=length_n,
                                     length_k=length_k, log=log)
            if not recursiv:
                break
    else:
        classical_matrix = csv_to_matrix(mypath, True)

        markovPatternAnalyse(mypath, classical_matrix, length_n=length_n, length_k=length_k, log=log)


def markovPatternAnalyse(path, classical_matrix, length_n, length_k, log):
    """
    from matrix to expression analysis
    :param path: path to file (necessary for the output name)
    :param classical_matrix: input matrix
    :param length_n: length of the tuple's to analyze
    :param length_k: length of the tuple's to extend
    :param log: if True use log scale
    :return:
    """
    shrinkSize = length_n
    classical_matrix_shrinked = matrix_shrink_size(shrinkSize, classical_matrix)
    if (length_k == 0):
        plot_image(classical_matrix_shrinked,
                   path[:-4] + "_k" + str(length_k) + "_markov-chain-approx_n" + str(shrinkSize),
                   sequence_size=np.sum(classical_matrix))
        representation_matrix = classical_matrix_shrinked
    else:
        representation_matrix = csv_representation_level_markov_chain(classical_matrix_shrinked,
                                                                      length_k=length_k)
        print(representation_matrix)
        print(classical_matrix)
        plot_image(representation_matrix,
                   path[:-4] + "_k" + str(length_k) + "_markov-chain-approx_n" + str(shrinkSize),
                   sequence_size=np.sum(classical_matrix))
    sorting(representation_matrix, shrinkSize,
            path[:-4] + "_k" + str(length_k) + "_n" + str(shrinkSize), False, log=log)


def path_to_fastaFiles(mypath, recursiv, n=8, filter_mikroSats=False):
    """
    parse multiple fasta files to csv
    :param mypath: path to fasta files
    :param recursiv: if True use also sub-dirs
    :param n: length of the tuple's
    :return:
    """
    global data
    if (os.path.isdir(mypath)):
        for root, dirs, files in os.walk(mypath):
            print(root)
            for filename in [f for f in files if
                             not f.endswith(".csv") and not f.endswith(".txt") and not f.endswith(".pdf")]:
                print(filename)
                address = str(os.path.join(root, filename))
                data = address
                matrix = readDNA(data, n, subseq=False, filter_mikroSats=filter_mikroSats)
                matrix_to_csv(matrix, address)
            if not recursiv:
                break
    elif (os.path.isfile(mypath)):
        data = mypath
        matrix = readDNA(data, n, subseq=False, filter_mikroSats=filter_mikroSats)
        matrix_to_csv(matrix, mypath)
    else:
        print("no file or directory!")
        print("change the path.")


def distance(seq, j, length):
    return editdistance.eval(seq[j: j + length], seq[j - length: j])


def find_microsats(dnaSeq, min_hits=10, max_len_tuple=5):
    """
    1. go over sequence
    2. try to extend pattern
    3. jump forward to next possible pattern
    :param dnaSeq:
    :return:
    """
    stp = 0
    microsats = []

    """1. go over sequence"""
    while stp < len(dnaSeq):
        ends = np.array([], dtype=int)
        hits = np.array([], dtype=int)
        length_sat = np.array([], dtype=int)

        """check all possible microsatellite tuple lengths"""
        for length in range(2, max_len_tuple + 1):
            pattern = dnaSeq[stp:stp + length]
            posi = stp + length
            mis = 0
            endposi = posi - 1
            hit = 1
            stp_2 = stp

            """2. try to extend pattern"""
            while posi < len(dnaSeq):
                dis = editdistance.eval(pattern, dnaSeq[posi: posi + length])

                # # check if substitution
                if dis >= 1 and dis < length:
                    dis_short = editdistance.eval(pattern, dnaSeq[posi: posi + length - 1])
                    if dis_short == 1:
                        dis_shift = editdistance.eval(pattern, dnaSeq[posi + length - 1:posi + length - 1 + length])
                        dis_next = editdistance.eval(pattern, dnaSeq[posi + length: posi + length + length])
                        if dis_shift < dis_next:
                            dis = 1
                            posi += -1

                if dis <= 1:
                    mis += dis
                    if mis <= 2:
                        hit += 1
                        endposi = posi + length - 1
                else:
                    # check for reading frame shift about 1
                    if (posi + length + 1 < len(dnaSeq)):
                        if (editdistance.eval(pattern, dnaSeq[posi + 1: posi + length + 1]) + mis <= 2):
                            # if insertion change reading frame
                            posi += 1 - length
                            mis += 1
                        else:
                            mis += dis

                    # check for reading frame shift about 2
                    elif (posi + length + 2 < len(dnaSeq) and mis == 0):
                        if (editdistance.eval(pattern, dnaSeq[posi + 2: posi + length + 2]) + mis <= 2):
                            # if indel change reading frame
                            posi += 2 - length
                            mis += 2
                        else:
                            mis += dis

                    else:
                        mis += dis
                posi += length

                if mis > 2:
                    # check if possible to extend to the left
                    if mis - dis <= 1 and stp - length + 1 >= 0:
                        dis_start = editdistance.eval(pattern, dnaSeq[stp - length:stp])
                        if dis_start + mis - dis <= 2:
                            stp_2 += -length
                            hit += 1
                        else:
                            dis_start = editdistance.eval(pattern, dnaSeq[stp - length + 1:stp])
                            if dis_start + mis - dis <= 2:
                                stp_2 += -length + 1
                                hit += 1
                    break

            ends = np.append(ends, endposi)
            hits = np.append(hits, hit)
            length_sat = np.append(length_sat, endposi - stp_2)

        if max(hits) >= min_hits:
            # find longest sequences of all sequences with max(hits)
            longestSeq = max(length_sat[hits >= min_hits])
            allowed_by_hits = hits >= min_hits
            allowed_by_length = length_sat == longestSeq
            # get end position of this sequence
            end = ends[allowed_by_hits * allowed_by_length]
            # assert len(end) == 1, str(allowed_by_hits * allowed_by_length)+"\n end:"+str(end)
            # if there are multiple entries then they have the same value
            end = int(end[0])

            # sat = (stp, max(ends))
            sat = (end - longestSeq, end)
            microsats.append(sat)
            """3. jump forward to next possible pattern"""
            # stp = max(ends)
            stp += 1
        else:
            stp += 1
    return microsats


def del_microsats(dnaSeq, min_hits=5, max_len_tuple=5):
    """
    1. get positions of microsatellites
    2. delete satellites
    :param dnaSeq: input sequence
    :param min_hits: minimal pattern hits to be microsat.
    :return: clean dna
    """

    """1. get positions of microsatellites"""
    microsats = find_microsats(dnaSeq, min_hits=min_hits, max_len_tuple=max_len_tuple)

    """2. reverse list for easier deletion"""
    microsats = list(reversed(microsats))

    """3. delete satellites"""
    dnaSeq = list(dnaSeq)

    for sat in microsats:
        # print(dnaSeq[sat[0]:sat[1] + 1])
        # del dnaSeq[sat[0]:sat[1] + 1]
        dnaSeq[sat[0]:sat[1] + 1] = map(str.lower, dnaSeq[sat[0]:sat[1] + 1])

    dnaSeq = ''.join(dnaSeq)
    return re.findall(r"[^acgt]+",dnaSeq)
    # return ''.join(dnaSeq)
    # return ''.join(x for x in dnaSeq if not x.islower())


# with open('/home/go96bix/Dropbox/hiwiManja/Fraktale/FractalDNA/test.fa') as f:
#     reg = re.compile(r">.*\n", re.IGNORECASE)
#     seq = f.readlines()
#     dnaSeq = str("".join(seq))
#     dnaSeq = reg.sub("", dnaSeq)
#     dnaSeq = dnaSeq.replace("\n", "").replace('\r', '')
#     dnaSeq = dnaSeq.upper()
#     dnaSeq = re.sub(r"[^ACGT]+", '', dnaSeq)
#     import time
#
#     start = time.process_time()
#     for i in range(1):
#         print(del_microsats(dnaSeq, min_hits=10))
#     end = time.process_time()
#     print(end - start)
# seq = "ccacacagatttt".upper()
# seq = "ccacacagacatcataacaaaaaattctcatttctatggtatgcccacacagacatcataacaaaaaattctcatttctatggtatgcccacacagacatcataacaaaaaattctcatttctatggtatgcccacacagacatcataacaaaaaattctcatttctatggtatgc".upper()
#
# print(del_microsats(seq,min_hits=5))
# readDNA("/home/go96bix/Dropbox/hiwiManja/Fraktale/FractalDNA/dinuShuffle.fa", 8, filter_mikroSats=True)
