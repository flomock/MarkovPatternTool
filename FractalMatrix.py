import re
# import sys
import math
import matplotlib

# Force matplotlib to not use any Xwindows backend.
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import os
from os import listdir
from os.path import isfile, join


c_score = 0
g_score = 0
t_score = 0
gc_score = 0
plt.rcParams['figure.figsize'] = 10, 10  # wi

def one_hot_encoding(seq):
    global g_score
    global c_score
    global t_score

    seq_array = np.array(list(seq))
    x_bin = np.zeros(seq_array.shape,dtype=np.int8)
    y_bin = np.zeros(seq_array.shape,dtype=np.int8)
    z_bin = np.zeros(seq_array.shape,dtype=np.int8)

    # x_bin[seq_array == 'A'] = 0
    # x_bin[seq_array == 'C'] = 1
    # c_score = x_bin.sum()
    # x_bin[seq_array == 'G'] = 1
    # g_score = x_bin.sum()-c_score
    # x_bin[seq_array == 'T'] = 0

    y_bin[seq_array == 'A'] = 1
    a_score = y_bin.sum()
    # y_bin[seq_array == 'C'] = 0
    y_bin[seq_array == 'G'] = 1
    # y_bin[seq_array == 'T'] = 0
    g_score = y_bin.sum() - a_score

    z_bin[seq_array == 'A'] = 1
    z_bin[seq_array == 'C'] = 1
    # z_bin[seq_array == 'G'] = 0
    # z_bin[seq_array == 'T'] = 0
    c_score = z_bin.sum() - a_score
    t_score = len(seq) - a_score - c_score - g_score

    return x_bin,y_bin,z_bin

def buildMatrix(size):
    # global matrix
    matrix = np.zeros((2 ** size, 2 ** size), dtype=np.int)
    return matrix


# einlesen der DNA und Berechung starten
def readDNA(startsize, endsize, subseq):
    global gc_score
    global c_score
    global g_score
    global t_score
    gc_score = 0
    matrixList = []
    for size in range(startsize, endsize + 1):
        matrixList.append(buildMatrix(size))

    if subseq == False:

        with open(data) as f:
            #    with open('shift4.mfa') as f:
            reg = re.compile(r">.*\n", re.IGNORECASE)
            seq = f.readlines()
            dnaSeq = str("".join(seq))
            dnaSeq = reg.sub("", dnaSeq)
            dnaSeq = dnaSeq.replace("\n", "").replace('\r', '')
            dnaSeq = dnaSeq.upper()
            dnaSeq = re.sub(r"[^ACGT]+", '', dnaSeq)
    else:
        dnaSeq = str(re.sub(r"[^ACGT]+", '', str(subseq)))


    x_binary, y_binary, z_binary = one_hot_encoding(dnaSeq)

    def bin_to_int(binary, bigest_bin, old = -1):
        """
        parse binary numbers to int and use fact that old already calculated
        :param binary:
        :param bigest_bin:
        :param old:
        :return:
        """
        decimal = old % bigest_bin
        decimal = decimal * 2 + int(binary[-1])
        if(old==-1):
            decimal = 0
            for digit in binary:
                decimal = decimal * 2 + int(digit)
        return decimal

    bigest_bin = 2**(endsize-1)
    old_y = -1
    old_z = -1

    for posi in range(0,len(dnaSeq),1):
        y = bin_to_int(y_binary[posi:posi + endsize],bigest_bin,old_y)
        z = bin_to_int(z_binary[posi:posi + endsize],bigest_bin,old_z)

        old_y = y
        old_z = z

        matrixList[0][y][z] += 1
    c_score = float(c_score) / float(len(dnaSeq))
    g_score = float(g_score) / float(len(dnaSeq))
    t_score = float(t_score) / float(len(dnaSeq))

    return matrixList


def identify(seq, j, length):
    if (j - length < 0):
        return False
    elif (seq[j: j + length] == seq[j - length: j]):
        return True
    else:
        return False


# ACHTUNG gilt nur fuer num i = y und num j = z
def int_to_binarray(num_i, num_j, size):
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

    ### To use if arrows wanted ###

    # cmap = plt.cm.gray_r
    ## norm = plt.Normalize(0,np.asarray(pic).max()) # vll eleganter als oberen 4 zeilen
    # rgba = cmap(pic)
    ## pixel_sequence = np.empty([1,2])get_pixel_sequence("/mnt/fass1/genomes/Eukaryots/homo_sapiens_done/all.hsa.ncbi.fa")
    # pixel_sequence = get_pixel_sequence("/home/go96bix/Dropbox/hiwiManja/Fraktale/FractalDNA/fasta/hs_ref_GRCh38.p2_chr21.mfa")

    # for k in range(0,len(pixel_sequence[:,0]),1):
    #     rgba[pixel_sequence[k,0],pixel_sequence[k,1],:3] = 1,0,0
    ## rgba[255,0,:3] = 1,0,0
    # plt.imshow(rgba, cmap='gray_r', interpolation="none", extent=[0, 1, 1, 0])


    ### Use if no arrwos wanted ###
    plt.imshow(pic, cmap='gray_r', interpolation="none", extent=[0, 1, 1, 0])
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

    # ax = show_path_of_tuble(ax,pixel_sequence) # if arrows wanted
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(outputName + ".pdf")
    plt.close()


def show_path_of_tuble(ax, pixel_sequence):
    for k in range(1, len(pixel_sequence[:, 0]), 1):
        # ,,:3] = 1,0,0
        former_x = float(pixel_sequence[k - 1, 1]) / 256.0
        former_y = float(pixel_sequence[k - 1, 0]) / 256.0

        ax.arrow(former_x, former_y,
                 (float(pixel_sequence[k, 1]) / 256.0) - former_x, (float(pixel_sequence[k, 0]) / 256.0) - former_y
                 , length_includes_head=True, head_width=0.01, head_length=0.02, fc='b', ec='r')
    return ax


def int_to_seq(num_i, num_j, size):
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


def sorting(outMatrix, size, outputName, ploting):
    # sorting
    sorted_seq = []
    m = outMatrix.flatten()
    m_sorted = np.argsort(m)
    if ploting == True:
        m_small = m_sorted[0:15]
        m_small = np.append(m_small, m_sorted[len(m_sorted) - 15:len(m_sorted)])
        m_sorted = m_small

    # calculate placement in array
    for number in m_sorted:
        column = number % (2 ** size)
        row = int((number - column) / (2 ** size))
        # seq = [(float(m[number])/float(np.sum(outMatrix))), int_to_seq(row, column, size)]
        seq = [m[number], int_to_seq(row, column, size)]
        sorted_seq.append(seq)

    if ploting == True:
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
        # plt.show()
    elif (ploting == "energy"):
        return sorted_seq
    else:
        # output file
        # print sorted_seq
        if "/" not in outputName:
            with open('/home/go96bix/Dropbox/hiwiManja/Fraktale/FractalDNA/MostLeast/' + outputName + ".txt",
                      "w") as text_file:
                for item in sorted_seq:
                    text_file.write("%s" % item[0] + "\t" + str("".join(item[1])) + "\n")
                    # np.savetxt(outputName,str(sorted_seq))
        else:
            with open(outputName + ".txt", "w") as text_file:
                for item in sorted_seq:
                    text_file.write("%s" % item[0] + "\t" + str("".join(item[1])) + "\n")


def matrix_to_csv(outMatrix, StandartOutputName):
    r = re.compile(r"/")
    outputName = r.sub("-", data) + ".csv"
    if StandartOutputName != True:
        # outputName = r.sub("-", data) + str(StandartOutputName) + ".csv"
        outputName = str(StandartOutputName) + ".csv"
    np.savetxt(outputName, outMatrix, fmt='%.1d', delimiter=',')


def csv_to_matrix(file, scoring):
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

def csvs_to_plot(mypath, size, pushContrast):
    paths = [mypath + '/' + f for f in listdir(mypath) if isfile(join(mypath, f))]
    for i in paths:
        foo = csv_to_matrix(i, True)
        name = str(i)[len(mypath) + 1:-4]
        plot_image(foo,name)


def csvs_to_moLeCommon(mypath, size, plotting):
    # paths = [mypath + '/' + f for f in listdir(mypath) if isfile(join(mypath, f))]
    # for i in paths:
    #     foo = csv_to_matrix(i, False)
    #     name = str(i)[len(mypath) + 1:-4]
    for root, dirs, files in os.walk(mypath):
        # print root
        for filename in [f for f in files if f.endswith(".csv")]:
            # print filename
            address = str(os.path.join(root, filename))
            foo = csv_to_matrix(address, False)
            if plotting:
                if not os.path.isdir(str(root + "/MostLeast/")):
                    os.makedirs(str(root + "/MostLeast/"))
                sorting(foo, size, str(root + "/MostLeast/" + filename[:-4]), plotting)
            else:
                sorting(foo, size, address[:-4], plotting)



def csv_representation_level(length_n, length_m, input_matrix):
    matrix_n = matrix_shrink_size(length_n, input_matrix)
    # plot_image(matrix_n,"test-n")
    if (length_n != length_m):
        matrix_m = matrix_shrink_size(length_m, input_matrix)
        # plot_image(matrix_m,"test-m")
    else:
        matrix_m = matrix_n
    prediction_matrix = buildMatrix(length_n + length_m).astype(float)
    # print matrix_n
    # print matrix_m
    total_n = np.sum(matrix_n)

    # total_m = np.sum(matrix_m)
    # print total_m
    # print total_n
    for i in range(0, 2 ** length_n):
        for j in range(0, 2 ** length_n):
            for k in range(0, 2 ** length_m):
                for l in range(0, 2 ** length_m):
                    prediction = float(matrix_n[i, j]) / float(total_n) * float(matrix_m[k, l]) / float(total_n)
                    # blaehe jedes quadrat in matrix n auf. Auf groesse [2**length_m,2**length_m]
                    # print prediction
                    prediction_matrix[(i * (2 ** length_m)) + k][(j * (2 ** length_m)) + l] = prediction
                    # print (str((i*(2**length_m))+k) +" , "+str((j*(2**length_m))+l))
    plot_image(prediction_matrix, "prediction")
    if (length_m + length_n != int(math.log(len(input_matrix[0]), 2))):
        out_matrix = matrix_shrink_size(length_m + length_n, input_matrix).astype(float) / total_n
        # out_matrix
    else:
        out_matrix = input_matrix.astype(float) / total_n
    # print out_matrix
    # print prediction_matrix
    out = [[float(out_matrix[row][col]) / float(prediction_matrix[row, col]) for row in range(len(out_matrix[0]))] for
           col in
           range(len(out_matrix[0]))]
    print(np.mean(out))

    # print out
    return np.asarray(out)
    # for i in range(0,2**(length_m+length_m)):
    #     for j in range(0,2**(length_m+length_m)):
    #         out_matrix


def csv_representation_level_markov_chain(input_matrix, length_k=2):
    """main method of pattern prediction with markov chains
    1. determine the probabilities for tuple appending
    e.g. how likely ACT --> ACTG
    2. make prediction how matrix with length n will look like
    3. return real- / predicted- data
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

    # helper_matrix = matrix_k[11,3]*prediction_matrix[(11-(2**(k-1)))*2][(3-(2**(length_k-1)))*2] #if 11 - 2**  >0 else 11 *2
    # dann aufblaehen siehe oben, aber nur um 2*2

    # for l in range(0, 2 ** size):
    #     for m in range(0, 2 ** size):
    # prediction_matrix[l,m]=matrix_k*transition_matrix #koennte rekurrenz (schwer)


def matrix_shrink_size(new_size, input_matrix):
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


### TEST ###

testcircle = -10
if testcircle == -1:
    # set path to fasta File in terminal like: python FractalMatrix.py /home/test.fa
    # if wanted to set file here in the code include next line and change the path
    # data "your/path/file.fa"
    matrixList = readDNA(8, 8, False)
    plot_image(matrixList[0], "test.csv")

if testcircle == 0:
    mypath = '/home/go96bix/Dropbox/hiwiManja/Fraktale/FractalDNA/treedata/fungi/fungi-markov'
    # presentation_2D(8,csv_to_matrix("/home/go96bix/Dropbox/hiwiManja/Fraktale/FractalDNA/csv types/precursor_RNA.csv",False),False,"precursor_RNA")

    # csvs_to_phyltree(mypath, 8, 1)
    # csvs_to_plot(mypath,8,False)
    # csvs_to_moLeCommon(mypath, 8, False)



# start program 2D Plot

if testcircle == 2:
    test_size = 8
    # dnaSeq = 0
    # with open("hs_ref_GRCh38.p2_chr21.mfa") as f:
    #     #    with open('shift4.mfa') as f:
    #     reg = re.compile(r">.*\n", re.IGNORECASE)
    #     seq = f.readlines()
    #     dnaSeq = bytearray("".join(seq))
    #     dnaSeq = reg.sub("", dnaSeq)
    #     dnaSeq = dnaSeq.replace("\n", "").replace('\r', '')
    #     dnaSeq = re.sub(r"[^ACGTacgt]+", '', dnaSeq)
    # matrix = readDNA(test_size, test_size, str(Seq(str(dnaSeq)).reverse_complement()))[0]
    # matrix_to_csv(matrix)
    # presentation_2D(test_size, matrix, False, data + "_reverse")
    # mypath = 'C:\Users\Florian\Dropbox\hiwiManja\Fraktale\FractalDNA\csvToPlot'
    mypath = '/home/go96bix/Dropbox/hiwiManja/Fraktale/FractalDNA/csvToPlot/human'
    for root, dirs, files in os.walk(mypath):
        print(root)
        for filename in [f for f in files if f.endswith(".csv")]:
            print(filename)

            if not os.path.isfile(str(os.path.join(root, str(filename[:-4] + ".pdf")))):
                var = input("Please enter testsize: ")
                test_size = int(var)
                if test_size == 0:
                    continue
                address = str(os.path.join(root, filename))
                plot_image(csv_to_matrix(address, True), address[:-4])

                # paths = [mypath + '/' + f for f in listdir(mypath) if isfile(join(mypath, f))]
                # for path in paths:
                #     presentation_2D(test_size, csv_to_matrix(path, True),False, str(data))
                #

                # start sorted plot
                # readDNA(8)
                # sorting(matrix,test_size)


# simply get csv of seq
if testcircle == 7:
    startsize = 8
    stopsize = 8
    # mypath = '/home/go96bix/Dropbox/hiwiManja/Fraktale/FractalDNA/gbk_fa.out'
    # paths = [mypath + '/' + f for f in listdir(mypath) if isfile(join(mypath, f))]
    with open("/home/go96bix/Dropbox/hiwiManja/Fraktale/FractalDNA/sequences/sequences")as f:
        # with open("C:\Users\Florian\Dropbox\hiwiManja\Fraktale\FractalDNA\sequences\sequencesRandom")as f:
        for line in f:
            data = line[:-1]  # ginge auch mit .strip['\n']
            if data[0] == '#':
                continue
                #  for file in paths:
                #      data=file
            k = -1 + startsize
            matrixlist = readDNA(startsize, stopsize, False)
            for m in matrixlist:
                k += 1
                matrix_to_csv(m, "cpg")




if testcircle == 13:
    size = 8
    mypath = '/home/go96bix/Dropbox/hiwiManja/Fraktale/FractalDNA/csvToPlot'
    for root, dirs, files in os.walk(mypath):
        print(root)
        for filename in [f for f in files if f.endswith(".csv")]:
            print(filename)

            if not False:  # os.path.isfile(str(os.path.join(root, str(filename[:-4] + ".pdf")))):
                var = input("Please enter testsize: ")
                test_size = int(var)
                if test_size == 0:
                    continue
                address = str(os.path.join(root, filename))
                classical_matrix = csv_to_matrix(address, True)
                prob_matrix = np.zeros((2 ** size, 2 ** size), dtype=np.float)

                # all posible combinations
                # get an array with length size and all combinations of 1,2,3
                a_score = 1 - (gc_score + t_score)

                for i in range(0, (4 ** size)):
                    wort = i
                    tupel = []
                    for j in range(0, size):
                        tupel.append(0)

                    for x in range(size - 1, -1, -1):
                        buchstabe = wort % 4
                        wort /= 4
                        wort = int(wort)
                        tupel[x] = buchstabe
                    # calcCoord(tupel, threeD_on, preview)
                    prob = 1.0
                    setOfLetters = []

                    for j in tupel:
                        # print j
                        if j == 0:
                            prob *= a_score
                        if j == 1:
                            prob *= c_score
                        if j == 2:
                            prob *= g_score
                        if j == 3:
                            prob *= t_score

                        setOfLetters.append(calcCoord(j))
                    # print prob
                    tupel = np.array(setOfLetters)
                    y = 0
                    z = 0
                    for k in range(0, size):
                        y += tupel[k, 1] * 2 ** (len(tupel[:, 1]) - (k + 1))
                        z += tupel[k, 2] * 2 ** (len(tupel[:, 2]) - (k + 1))
                    prob_matrix[y][z] = prob
                # print prob_matrix
                print(
                    "g_score " + str(g_score) + "c_score " + str(c_score) + "t_score " + str(
                        t_score) + "a_score: " + str(
                        a_score))
                result_matrix = np.divide(classical_matrix, prob_matrix)
                # decide which plots/csv zou need.
                # probably the last plot_image is what zou need

                plot_image(classical_matrix, address[:-4])
                # plot_image(prob_matrix, address[:-4] + "-probability")
                # presentation_2D(test_size, classical_matrix, True, address[:-4])
                # presentation_2D(8,prob_matrix,True,address[:-4]+"-probability")
                matrix_to_csv(result_matrix, address[:-4] + "-withProb")
                plot_image(result_matrix, address[:-4] + "-withProb",
                           sequence_size=int((round(np.sum(classical_matrix) / 100)) * 100))

def path_to_markovPatternAnalyse(mypath,length_n,length_k,recursiv):
    if (os.path.isdir(mypath)):
    # mypath = '/home/go96bix/Dropbox/hiwiManja/Fraktale/FractalDNA/csvToPlot'
        for root, dirs, files in os.walk(mypath):
            # print root
            for filename in [f for f in files if f.endswith(".csv")]:
                # print filename
                address = str(os.path.join(root, filename))
                classical_matrix = csv_to_matrix(address, True)
                markovPatternAnalyse(path=address, classical_matrix=classical_matrix,length_n=length_n,length_k=length_k)
            if not recursiv:
                break
    else:
        classical_matrix = csv_to_matrix(mypath, True)

        markovPatternAnalyse(mypath,classical_matrix,length_n=length_n,length_k=length_k)

def markovPatternAnalyse(path, classical_matrix, length_n, length_k):
    shrinkSize = length_n
    classical_matrix_shrinked = matrix_shrink_size(shrinkSize, classical_matrix)
    if (length_k == 0):
        plot_image(classical_matrix_shrinked, path[:-4] + "_k" + str(length_k) + "_markov-chain-approx_n" + str(shrinkSize),
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
    # matrix_to_csv(representation_matrix,address[:-4]+"_k"+str(length_k)+"_markov-chain-approx")
    ## plot_image(representation_matrix,address[:-4] + "-NEW")
    sorting(representation_matrix, shrinkSize,
            path[:-4] + "_k" + str(length_k) + "_n" + str(shrinkSize), False)

def path_to_fastaFiles(mypath,recursiv,n=8):
    global data
    if (os.path.isdir(mypath)):
        for root, dirs, files in os.walk(mypath):
            print(root)
            for filename in [f for f in files if not f.endswith(".csv") and not f.endswith(".txt") and not f.endswith(".pdf")]:
                print(filename)
                address = str(os.path.join(root, filename))
                data = address
                matrixlist = readDNA(n,n,False)
                for m in matrixlist:
                    matrix_to_csv(m,address)
            if not recursiv:
                break
    elif(os.path.isfile(mypath)):
        data = mypath
        matrixlist = readDNA(n, n, False)
        for m in matrixlist:
            matrix_to_csv(m, mypath)
    else:
        print("no file or directory!")
        print("change the path.")
