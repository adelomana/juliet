import sys, os, time

def memeReader(module, outputFileForMeme):

    # define some variables
    motifs = []
    transcriptHits = []

    # read the information
    reading = False
    with open(outputFileForMeme, 'r') as f:
        for line in f:
            vector = line.split()
            if len(vector) == 7:
                if vector[6] == 'p-value':
                    reading = True
                if reading == True and vector[0] == 'Motif':
                    motifs.append([]); transcriptHits.append([])
                if reading == True and vector[0][:3] == 'Cre':
                    pvalue = float(vector[3])
                    sequence = vector[5]
                    transName = vector[0]
                    if motifs[-1] == []:
                        motifs[-1] = [sequence]; transcriptHits[-1] = [transName]
                    else:
                        if sequence == motifs[-1][0] and transName != transcriptHits[-1][0]:
                            motifs[-1].append(sequence); transcriptHits[-1].append(transName)

    #print motifs
    #print transcriptHits

    return motifs, transcriptHits

def complementarySequenceConverter(a):

    if a == '':
        print('error 2 from complementarySequenceConverter')
        sys.exit() 

    b = ''
    for x in a:
        if x == 'A':
            c = 'T'
        elif x == 'T':
            c = 'A'
        elif x == 'G':
            c = 'C'
        elif x == 'C':
            c = 'G'
        elif x == 'N':
            c = 'N'
        else:
            print('error 1 from complementarySequenceConverter')
            print(c)
            sys.exit()
        b = b + c

    return b

def memeRunner(module, transcripts, sequencesSet, numberOfMotifs, motifLength, memeDir):
    
    inputFileForMeme = '%s/inputFileForMEME_%s.txt'%(memeDir, module)
    outputFileForMeme = '%s/resultsMEME_%s.txt'%(memeDir, module)
    errorFileForMeme = '%s/errorMessagesFromMEME.txt'%memeDir
    
    # create input for meme
    g=open(inputFileForMeme, 'w')
    for i in range(len(transcripts)):
        g.write('>%s\n'%transcripts[i])
        g.write('%s\n'%sequencesSet[i])
    g.close()

    # execute meme
    cmd='time meme %s -dna -revcomp -mod anr -w %s -nmotifs %s -text > %s 2> %s'%(inputFileForMeme, motifLength, numberOfMotifs, outputFileForMeme, errorFileForMeme)
    print(cmd)
    os.system(cmd)

    return outputFileForMeme

def reader():

    sourceDir='/Volumes/~alomana/projects/green/results/cMonkey/boyle12/run.36/CMONKE~1.21_'
    GRMs={}
    inputFile=sourceDir+'/cluster.members.genes.txt'
    with open(inputFile) as f:
        for line in f:
            vector=line.split()
            localName=int(vector[0])
            localTranscripts=vector[1:len(vector)]
            GRMs[localName]=localTranscripts

    # printing low confidence network
    lowConfidenceTranscripts=[]
    for module in GRMs:
        for element in GRMs[module]:
            if element not in lowConfidenceTranscripts:
                lowConfidenceTranscripts.append(element)
    print('low confidence transcripts',len(lowConfidenceTranscripts))
    print('low confidence modules',len(GRMs))

    # cleaning out the ones that have residual higher than 0.4 and no motif below 0.1
    GRMsQuality={}
    inputFile=sourceDir+'/cluster.summary.tsv'
    with open(inputFile) as f:
        next(f)
        for line in f:
            vector=line.split('\t')
            name=int(vector[0])
            residual=float(vector[4])
            eval1=float(vector[6])
            eval2=float(vector[8])
            eval3=min([eval1,eval2])
            if residual > 0.4 and eval3 > 0.1:
                del GRMs[name]

    # printing high confidence network
    highConfidenceTranscripts=[]
    for module in GRMs:
        for element in GRMs[module]:
            if element not in highConfidenceTranscripts:
                highConfidenceTranscripts.append(element)
    print('high confidence transcripts',len(highConfidenceTranscripts))
    print('high confidence modules',len(GRMs))
                                
    return GRMs

def sequencesRetriever(transcripts, upstreamLength, downstreamLength):

    sequencesSet = []

    for transcript in transcripts:
        upstream = upstreamRetriever(transcript, upstreamLength, downstreamLength)
        sequencesSet.append(upstream)

    return sequencesSet

def upstreamRetriever(transcript, upstreamLength, downstreamLength):

    # defining the name of the gene, the position in the genome and the strand 
    geneName = transcript.split('.t')[0]
    found = False

    gff3File='/Users/alomana/projects/green/data/genomic/phytozome/v5.5/annotation/Creinhardtii_281_v5.5.gene.gff3'
    with open(gff3File, 'r') as f:
        next(f)
        next(f)
        for line in f:
            vector = line.split('\t')
            if vector[2] == 'gene':
                possible=vector[8].split('.v5.5')[0].split('ID=')[1]
                if possible == geneName:
                    chrom = vector[0]
                    start = int(vector[3])
                    end = int(vector[4])
                    strand = vector[6]
                    found = True

    # retrieving the chromosome sequence
    reading = False
    chromSeq = ''
    genomeFile = '/Users/alomana/projects/green/data/genomic/phytozome/v5.5/assembly/Creinhardtii_281_v5.0.fa'
    with open(genomeFile, 'r') as f:
        for line in f:
            line=line.replace('\n', '')

            if line[0] == '>':
                currentChromosome = line.split('>')[1]

            if chrom == currentChromosome:
                reading = True

            if reading == True and chrom != currentChromosome:
                break

            if reading == True and line[0] != '>':
                chromSeq = chromSeq + line

    # selecting the upstream depending on the strand, including a region downstream the transcriptional start site
    if strand == '+':
        tomos=(start - upstreamLength - 1, start - 1 + downstreamLength)
        upstream = chromSeq[tomos[0] : tomos[1]]
    elif strand == '-':
        tomos=(end - downstreamLength, end + upstreamLength)
        complementarySequence = chromSeq[tomos[0] : tomos[1]] 
        inverse = complementarySequence[::-1]
        upstream = complementarySequenceConverter(inverse)
        
    else:
        print('error 1 from upstreamRetriever')
        sys.exit()

    return upstream

### MAIN

print()
print('here comes the sun...')
print()

# 0. general variables

# 0.1 user defined variables
numberOfMotifs = 4 # should be 4 for real, 2 for testing
motifLength = 20
upstreamLength = 2000 # should be 2000 for real, 500 for testing
downstreamLength = 25

# 0.2 some paths
timestamp = '.'.join([str(element) for element in time.localtime()[0:6]])
#timestamp = 'test'
memeDir = 'meme_%s'%timestamp
if os.path.exists(memeDir) == False:
    os.mkdir(memeDir)
resultsFile = '%s/resultsFile.txt'%memeDir
res = open(resultsFile, 'w')
res.write('module\tnumberOfGenes\tmoduleSize\tmotif\n')

# 1. read the input of the structure
print('reading cMonkey output...')
dataStructure = reader()

# 2. for each module...
print()
for module in dataStructure:
    print('working on module', module)
    
    # 2.1 define the sequences
    print('\t retrieving upstream sequences...')
    transcripts = dataStructure[module]
    sequencesSet = sequencesRetriever(transcripts, upstreamLength, downstreamLength)

    # 2.2. run meme
    print('\t searching for common motifs...')
    outputFileForMeme = memeRunner(module, transcripts, sequencesSet, numberOfMotifs, motifLength, memeDir)
    
    motifs, transcriptHits = memeReader(module, outputFileForMeme)
    print(outputFileForMeme)
    print(motifs)
    print(transcriptHits)

    # 2.3. report the number of genes and the consensus sequence
    print('\t reading MEME results...')
    numberTargetGenes = [len(element) for element in motifs]
    sequences = [element[0] for element in motifs]

    for i in range(len(motifs)):
        if numberTargetGenes[i] >= 2:
            listOfHits = ','.join(transcriptHits[i])
            res.write('%s\t%s\t%s\t%s\t%s\n'%(module, numberTargetGenes[i], str(len(transcripts)), motifs[i][0], listOfHits))
            res.flush()
    print()
    
res.close()

print('... the world is a much better now.')
