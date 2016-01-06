import sys
import os
import multiprocessing
import math

# Path to BLAST
BLAST = "/blast-2.2.26/bin/"
# Path to nr
NR = "/nr/"
# Path to DISOPRED
DISOPRED = "/DISOPRED/"
# Path to Weka
WEKA = "/weka-3-6-11/"

def main():
    '''Runs the pipeline.'''
    # Get command-line parameters
    directory = 0
    threads = 1
    k = 0
    alpha = 0
    runblast = 1
    feats = "10010011101111"
    t = .5
    for i in range(len(sys.argv)):
        if (sys.argv[i] == "-d" and (i+1) < len(sys.argv)):
            directory = sys.argv[i+1]
        if (sys.argv[i] == "-n" and (i+1) < len(sys.argv)):
            threads = int(sys.argv[i+1])
        if (sys.argv[i] == "-k" and (i+1) < len(sys.argv)):
            k = int(sys.argv[i+1])
        if (sys.argv[i] == "-a" and (i+1) < len(sys.argv)):
            alpha = float(sys.argv[i+1])
        if (sys.argv[i] == "-b" and (i+1) < len(sys.argv)):
            runblast = int(sys.argv[i+1])
        if (sys.argv[i] == "-f" and (i+1) < len(sys.argv)):
            feats = int(sys.argv[i+1])
        if (sys.argv[i] == "-t" and (i+1) < len(sys.argv)):
            t = float(sys.argv[i+1])
    if (directory == 0):
        print("error: invalid parameters\n-d directory\n-n number of threads\n-k k neighbors for KNN\n-a alpha value for KNN\n-b run blast\n-f feature selection (see readme)\n-t threshold value")
        return 0
    if (k > 200):
        print("error: k cannot be greater then 200.")
        return 0
    # Run pipeline
    os.system("cd Data;rm -r -f " + directory + "_dis;rm -f " + directory + ".arff " + directory + "_binding.txt " + directory + "_data.txt " + directory + "_entropy.txt " + directory + "_props.txt trainer.arff predictions.txt;cd ../" + directory + ";ls | grep .fa > proteinList.txt")
    if (runblast):
        os.system("rm -r -f PSSM_" + directory)
        print("generating pssm matrices")
        generatePSSM(directory, threads)
    print("generating data file")
    generateData(directory)
    print("running KNN")
    KNN(directory, k, alpha, threads)
    print("extracting entropy")
    extractEntropy(directory)
    print("extracting properties")
    extractProperties(directory)
    print("extracting disorder")
    extractDisorder(directory)
    print("generating arff file")
    generateARFF(directory)
    print("generating trainer")
    generateTrainer(k, alpha)
    print("running Weka")
    runWeka(directory, feats)
    print("writing predictions")
    writePredictions(directory, t)
    print("done.")

def generatePSSM(directory, threads):
    '''Takes a directory name followed by the number of threads you want to use and runs PSI-blast on the queries listed in proteinList.txt against the nr-database with an e value of 0.001 for three iterations.'''
    # Initialize
    os.system("rm -r -f PSSM_" + directory + ";mkdir PSSM_" + directory + ";cd "+ directory)
    tempfile = open(directory + '/proteinList.txt', 'r')
    inventory = tempfile.readlines()
    tempfile.close()
    # Run BLAST and record results
    for query in inventory:
        os.system("cd " + directory + ';' + BLAST + "blastpgp -a " + str(threads) + " -d " + NR + "nr -i " + query[:-1] + " -j 3 -Q ../PSSM_" + directory + '/' + query.split('.')[0] + ".txt -h 0.001 > temp")
    matrixList = open("./PSSM_" + directory + "/matrixList.txt", 'w')
    for query in inventory:
        matrixList.write(query.split('.')[0] + ".txt\n")
    os.system("cd " + directory + ";rm -f error.log temp")
    matrixList.close()

def generateData(directory):
    '''Takes the name of a directory (and therefore also know the PSSM directory name) and makes a data file.'''
    # Intialize
    sequenceDictionary = {}
    tempfile = open(directory + "/proteinList.txt", 'r')
    proteinInventory = tempfile.readlines()
    tempfile.close()
    tempfile = open("PSSM_" + directory + "/matrixList.txt", 'r')
    pssmInventory = tempfile.readlines()
    tempfile.close()
    # Store data
    for pssm in pssmInventory:
        temp = open("PSSM_" + directory + '/' + pssm[:-1], 'r').readlines()
        modifiedMatrix = []
        for line in range(3, len(temp) - 6):
            modifiedMatrix.append(temp[line].split()[2:22])
            sequenceDictionary[pssm[:-5]] = (modifiedMatrix)
    # Manipulate and write data
    sequenceList = list(sequenceDictionary.keys())
    outputFile = open("Data/" + directory + "_data.txt", 'w')
    toWrite = ''
    for sequence in sequenceList:
        sequenceLength = len(sequenceDictionary[sequence]) - 1
        for protein in range(sequenceLength):
            toWrite += str(sequence) + '_' + str(protein) + " | "
            for i in range(-7, 8):
                index = protein + i
                if (index < 0):
                    index += sequenceLength
                elif (index >= sequenceLength):
                    index -= sequenceLength
                for item in sequenceDictionary[sequence][index]:
                    toWrite += item + ' '
                toWrite += "| "
            toWrite = toWrite[:-3] + '\n'
    outputFile.write(toWrite)    
    outputFile.close()

def KNN(directory, k, alpha, numThreads):
    '''Takes a value k and an alpha value and calculates the binding percentage for novel data in the specified directory.'''
    # Initialize
    temp = open("Data/" + directory + "_data.txt", 'r')
    test = temp.readlines()
    temp.close()
    temp = open("Data/" + "trainer.txt", 'r')
    trainer = temp.readlines()
    temp.close()
    results = open("Data/" + directory + "_binding.txt", 'w')
    # Save and interpret input
    trainingSet = []
    testSet = []
    for protein in trainer:
        trainingSet.append([protein.split('|')[0],[],protein.split('|')[16]])
        for profile in protein.split('|')[1:16]:
            trainingSet[-1][1].append(profile.split())
    for protein in test:
        testSet.append([protein.split('|')[0],[]])
        for profile in protein.split('|')[1:16]:
            testSet[-1][1].append(profile.split())
    # Para-Proteins
    sem = multiprocessing.Semaphore(numThreads)
    lock = multiprocessing.Lock()
    threadList = []
    for currentProtein in testSet:
        sem.acquire()
        lock.acquire()
        threadList.append(multiprocessing.Process(target=calc, args=(currentProtein, trainingSet, k, alpha, results, lock, sem)))
        threadList[-1].start()
        lock.release()
    for thread in threadList:
        thread.join()
    results.close()

def calc(currentProtein, trainingSet, k, alpha, results, lock, sem):
    '''Calculates the binding via KNN'''
    # Intialize variable for calculating distance
    distances = []
    neighbors = []
    # Loop through the training proteins
    for featureProtein in trainingSet:
        # Calculate the distance (sumation) and add it to a list of distances
        distance = 0
        for i in range(15):
            for j in range(20):
                distance += (8 - abs(7 - i)) ** 2 * pow(int(currentProtein[1][i][j]) - int(featureProtein[1][i][j]), 2)
        distances.append(math.sqrt(distance))
        # Is the neighbor a nearest neighbor?
        if (k > len(neighbors) or distances[-1] < neighbors[-1][0]):
            index = 0
            while (len(neighbors) > abs(index) and distances[-1] < neighbors[index - 1][0]):
                index -= 1
            if (index == 0):
                neighbors.append((distances[-1],featureProtein[0][:-1],featureProtein[2][-2]))
            else:
                neighbors.insert(index,(distances[-1],featureProtein[0][:-1],featureProtein[2][-2]))
            if (k < len(neighbors)):
                neighbors.pop()
    # Find average distance
    averageDistance = sum(distances) / len(distances)
    # Find standard deviation
    standardDeviation = 0
    for distance in distances:
        standardDeviation += (distance - averageDistance) ** 2
    standardDeviation = (standardDeviation / (len(distances) - 1)) ** (0.5)
    # Find zScores and write raw data to a results file
    toWrite = currentProtein[0] + " | "
    zscores = []
    for neighbor in range(len(neighbors)):
        zscores.append(((averageDistance - neighbors[neighbor][0])/standardDeviation, neighbors[neighbor][-1]))
    positiveBinding = 0
    negativeBinding = 0
    for zscore in zscores:
        if (zscore[1] == '1'):
            positiveBinding += float(zscore[0]) ** alpha
        else:
            negativeBinding += float(zscore[0]) ** alpha
    # Write output
    lock.acquire()
    results.write(currentProtein[0] + '\t' + str(positiveBinding / (positiveBinding + negativeBinding)) + '\n')
    results.flush()
    lock.release()
    sem.release()

def extractEntropy(data):
    '''Extracts entropy for secondary classification'''
    temp = open("PSSM_" + data + "/matrixList.txt", 'r')
    mlist = temp.readlines()
    temp.close()
    results = open("Data/" + data + "_entropy.txt", 'w')
    for m in mlist:
        temp = open("PSSM_" + data + '/' + m[:-1], 'r')
        current = temp.readlines()
        temp.close()
        for i in range(3,len(current)-6):
            e = 0
            for j in range(22,42):
                val = int(current[i].split()[j])
                if val > 0:
                    val /= 100
                    e += val * math.log(val, 20)
            results.write(m[:-5] + '_' + str(i-3) + '\t' + str(e) + '\n')
    results.close()

def extractProperties(data):
    '''Extracts atomic properties for secondary classification'''
    temp = open("Data/aaupdated", 'r')
    aai = temp.readlines()
    temp.close()
    temp = open(data + "/proteinList.txt", 'r')
    plist = temp.readlines()
    temp.close()
    results = open("Data/" + data + "_props.txt", 'w')
    for p in plist:
        temp = open(data + '/' + p[:-1], 'r')
        current = temp.readlines()
        temp.close()
        for i in range(len(current[1])-1):
            results.write(p[:-4] + '_' + str(i) + '\t' + str(i) + '\t' + current[1][i])
            if current[1][i] == 'U':
                for j in range(9):
                    results.write("\t0")
            else:
                populate(aai, current[1][i], results)
            results.write('\n')
    results.close()

def populate(aai, char, results):
    '''Helper to extract properties'''
    for i in range(len(aai)):
        if(aai[i][0] == 'H'):
            i += 1
            tlist = aai[i].split()[1:]
            for j in range(len(tlist)):
                if(char in tlist[j]):
                    if(char == tlist[j].split('/')[0]):
                        i += 1
                    else:
                        i += 2
                    results.write('\t' + str(aai[i].split()[j]))
                    results.flush()

def extractDisorder(data):
    '''Extracts disorder for secondary classification'''
    os.system("cp -r " + data + " Data/" + data + "_dis")
    temp = open("Data/" + data + "_dis/proteinList.txt", 'r')
    plist = temp.readlines()
    temp.close()
    for p in plist:
        os.system(DISOPRED + "run_disopred.pl Data/" + data + "_dis/" + p[:-1] + "> temp")
    os.system("rm -f temp")

def generateARFF(d):
    '''Generate the arff in order to run Weka'''
    p = open("Data/"+d+"_props.txt", "r").readlines()
    e = open("Data/"+d+"_entropy.txt", "r").readlines()
    aarf = open("Data/"+d+".arff", "w")
    aarf.write("@RELATION ATPbinding\n\n@ATTRIBUTE KNNprediction NUMERIC\n@ATTRIBUTE Position NUMERIC\n@ATTRIBUTE Residue {G,A,V,L,I,S,T,C,M,P,H,R,N,Q,E,D,F,W,Y,K,U}\n@ATTRIBUTE BULH740101 NUMERIC\n@ATTRIBUTE EISD840101 NUMERIC\n@ATTRIBUTE HOPT810101 NUMERIC\n@ATTRIBUTE RADA880108 NUMERIC\n@ATTRIBUTE ZIMJ680104 NUMERIC\n@ATTRIBUTE MCMT640101 NUMERIC\n@ATTRIBUTE BHAR880101 NUMERIC\n@ATTRIBUTE CHOC750101 NUMERIC\n@ATTRIBUTE COSI940101 NUMERIC\n@ATTRIBUTE Entropy NUMERIC\n@ATTRIBUTE Disorder NUMERIC\n@ATTRIBUTE Binding {0,1}\n\n@DATA\n")
    temp = open("Data/"+d+"_binding.txt", "r")
    b = temp.readlines()
    temp.close()
    p.sort()
    e.sort()
    b.sort()
    prev = ("",0,0,0)
    for c in b:
        c = c.split()
        if c[0][:4] == prev[0]:
            j = prev[1]
            k = prev[2]
        else:
            j = 0
            k = 0
        aarf.write(c[1]+",")
        while p[j] == "\n" or c[0] != p[j].split()[0]:
            j+=1
        t = p[j].split()
        aarf.write(t[1]+","+t[2]+","+t[3]+","+t[4]+","+t[5]+","+t[6]+","+t[7]+","+t[8]+","+t[9]+","+t[10]+","+t[11]+",")
        while c[0] != e[k].split()[0]:
            k+=1
        aarf.write(e[k].split()[1]+",")
        temp = open("Data/"+d+"_dis/"+c[0][:6]+".diso", "r")
        td = temp.readlines()
        temp.close()
        l = 0
        while len(td[l].split()) == 0 or td[l].split()[0] != str(int(c[0].split("_")[-1])+1):
            l+=1
        aarf.write(td[l].split()[-1]+",?\n")
        prev = (c[0][:4],j,k)

def generateTrainer(k, alpha):
    '''Generate the trainer arff in order to run Weka'''
    # Calculate binding
    b = []
    proteins = open("Data/trainer_zscores.txt", "r")
    for protein in proteins:
        temp = protein.split('|')
        proteinName = temp[0][:-2]
        binding = temp[-1][:-1]
        zscores = []
        for zscore in temp[1:-1]:
            zscores.append(zscore.split())
        positiveBinding = 0
        negativeBinding = 0
        for zscore in zscores:
            if (zscore[1] == '1'):
                positiveBinding += float(zscore[0]) ** alpha
            else:
                negativeBinding += float(zscore[0]) ** alpha
        b.append(proteinName + "\t" + str(positiveBinding / (positiveBinding + negativeBinding)) + "\t" + binding)
    # Write ARFF file
    p = open("Data/trainer_props.txt", "r").readlines()
    e = open("Data/trainer_entropy.txt", "r").readlines()
    aarf = open("Data/trainer.arff", "w")
    aarf.write("@RELATION ATPbinding\n\n@ATTRIBUTE KNNprediction NUMERIC\n@ATTRIBUTE Position NUMERIC\n@ATTRIBUTE Residue {G,A,V,L,I,S,T,C,M,P,H,R,N,Q,E,D,F,W,Y,K,U}\n@ATTRIBUTE BULH740101 NUMERIC\n@ATTRIBUTE EISD840101 NUMERIC\n@ATTRIBUTE HOPT810101 NUMERIC\n@ATTRIBUTE RADA880108 NUMERIC\n@ATTRIBUTE ZIMJ680104 NUMERIC\n@ATTRIBUTE MCMT640101 NUMERIC\n@ATTRIBUTE BHAR880101 NUMERIC\n@ATTRIBUTE CHOC750101 NUMERIC\n@ATTRIBUTE COSI940101 NUMERIC\n@ATTRIBUTE Entropy NUMERIC\n@ATTRIBUTE Disorder NUMERIC\n@ATTRIBUTE Binding {0,1}\n\n@DATA\n")
    p.sort()
    e.sort()
    b.sort()
    prev = ("",0,0,0)
    for c in b:
        c = c.split()
        if c[0][:4] == prev[0]:
            j = prev[1]
            k = prev[2]
        else:
            j = 0
            k = 0
        aarf.write(c[1]+",")
        while p[j] == "\n" or c[0] != p[j].split()[0]:
            j+=1
        t = p[j].split()
        aarf.write(t[1]+","+t[2]+","+t[3]+","+t[4]+","+t[5]+","+t[6]+","+t[7]+","+t[8]+","+t[9]+","+t[10]+","+t[11]+",")
        while c[0] != e[k].split()[0]:
            k+=1
        aarf.write(e[k].split()[1]+",")
        temp = open("Data/trainer_dis/"+c[0][:6]+".diso", "r")
        td = temp.readlines()
        temp.close()
        l = 0
        while len(td[l].split()) == 0 or td[l].split()[0] != str(int(c[0].split("_")[-1])+1):
            l+=1
        aarf.write(td[l].split()[-1]+","+c[-1]+"\n")
        prev = (c[0][:4],j,k)

def runWeka(d, feats):
    '''Runs Weka'''
    # Initialize variables
    predictions = ""
    roc = 0
    tn = 0
    fn = 0
    fp = 0
    tp = 0
    # Feature selection
    train = open("train"+feats+".arff","w")
    tempfile = open("Data/trainer.arff","r")
    temp = tempfile.readlines()
    tempfile.close()
    train.write(temp[0]+temp[1])
    for j in range(14):
        if feats[j] == "1":
            train.write(temp[j+2])
    train.write(temp[16]+temp[17]+temp[18])
    for j in range(19, len(temp)):
        line = temp[j].split(",")
        for k in range(14):
            if feats[k] == "1":
                train.write(line[k]+",")
        train.write(line[-1])
    train.close()
    test = open("test"+feats+".arff","w")
    tempfile = open("Data/"+d+".arff","r")
    temp = tempfile.readlines()
    tempfile.close()
    test.write(temp[0]+temp[1])
    for j in range(14):
        if feats[j] == "1":
            test.write(temp[j+2])
    test.write(temp[16]+temp[17]+temp[18])
    for j in range(19, len(temp)):
        line = temp[j].split(",")
        for k in range(14):
            if feats[k] == "1":
                test.write(line[k]+",")
        test.write(line[-1])
    test.close()
    # Run Weka
    os.system("export CLASSPATH="+WEKA+"weka.jar;java weka.classifiers.meta.ThresholdSelector -t train"+feats+".arff -T test"+feats+".arff -i -p 0 > temp"+feats)
    tempfile = open("temp"+feats,"r")
    temp = tempfile.readlines()[5:]
    tempfile.close()
    os.system("rm train"+feats+".arff test"+feats+".arff temp"+feats)
    for line in temp:
        line = line.split()
        if len(line):
            if line[2].split(":")[0] == "1":
                predictions += str(1 - float(line[-1])) + "\n"
            else:
                predictions += line[-1] + "\n"
    results = open("Data/predictions.txt", "w")
    results.write(predictions)
    results.close()
    
def writePredictions(d, t):
    '''Outputs the final predictions'''
    temp = open("Data/"+d+"_binding.txt", "r")
    b = temp.readlines()
    temp.close()
    b.sort()
    temp = open("Data/predictions.txt", "r")
    p = temp.readlines()
    temp.close()
    proteins = {}
    for line in b:
        c = line.split()[0]
        if c.split('_')[0]+'_'+c.split('_')[1] not in proteins.keys():
            proteins[c.split('_')[0]+'_'+c.split('_')[1]] = int(c.split('_')[-1])
        elif int(c.split('_')[-1]) > proteins[c.split('_')[0]+'_'+c.split('_')[1]]:
            proteins[c.split('_')[0]+'_'+c.split('_')[1]] = int(c.split('_')[-1])
    results = {}
    os.system('rm -r -f Results; mkdir Results')
    for key in list(proteins.keys()):
        results[key] = [0] * (proteins[key] + 1)
    for i in range(len(b)):
        c = b[i].split()[0]
        if float(p[i].split()[0]) >= t:
            results[c.split('_')[0]+'_'+c.split('_')[1]][int(c.split('_')[-1])] = "1"
        else:
            results[c.split('_')[0]+'_'+c.split('_')[1]][int(c.split('_')[-1])] = "0"
    for key in list(results.keys()):
        temp = open(d+"/"+key+".fa", "r")
        tf = temp.readlines()
        temp.close()
        cr = open("Results/"+key+".fa", "w")
        for line in tf[:-1]:
            cr.write(line)
        for binding in results[key]:
            cr.write(str(binding))
        cr.close()

main()
