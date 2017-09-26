
from sys import *

def parseInput(arg):

    if len(arg) != 4:
         print "wrong number of arguments; correct usage: python ./WP.py tree substitution_matrix leaf_seqs"
         exit(1)
    ###########################################################################
    try:
        tree = {}
        infile = open(arg[1], "r")
        line = infile.readline()
        line = line.strip()
        tree[line] = []
        root = line
        for entry in infile:
            entry = entry.strip()
            entry = entry.split(" ")
            if entry[0] not in tree:
                tree[entry[0]] = []
            if entry[1] in tree:
                tree[entry[1]].append(entry[0])
            else:
                tree[entry[1]] = [entry[0]]
            #tree[entry[0]] = []
    except:
        print "failed to read tree"
        exit(1)
    ###########################################################################
    try:
        subMatrix = {}
        infile = open(arg[2], "r")
        for line in infile:
            line = line.strip()
            items = line.split(" ")
            if items[0] in subMatrix:
                vals = subMatrix.get(items[0])
                if items[1] not in vals:
                    vals[items[1]] = float(items[2])
                    subMatrix[items[0]] = vals
            else:
                subMatrix[items[0]] = {}
                subMatrix[items[0]][items[0]] = 0.0
                subMatrix[items[0]][items[1]] = float(items[2])
            if items[1] in subMatrix:
                vals = subMatrix.get(items[1])
                if items[0] not in vals:
                    vals[items[0]] = float(items[2])
                    subMatrix[items[1]] = vals
            else:
                subMatrix[items[1]] = {}
                subMatrix[items[1]][items[1]] = 0.0
                subMatrix[items[1]][items[0]] = float(items[2])
    except:
        print "failed to read in substitution matrix"
        exit(1)
    ###########################################################################
    try:
        seqs= {}
        infile = open(arg[3], "r")
        for line in infile:
            line = line.strip()
            items = line.split(": ")
            seqs[items[0]] = items[1]
            seqlen = len(items[1])
    except:
        print "failed to read in sequence file"
    ###########################################################################

    return [tree, subMatrix, seqs, seqlen, root]

def wp(tree, subMatrix, seqs):
    #so for each position, we need to keep track of the scores for each possible option, and the children scores that lead to it
    #let's make a dictionary of dictionaries that does!

    scores = {}

    for node in tree:
        if len(tree[node]) == 0:
            #finding the leaf nodes to set their values
            scores[node] = leafNode(node, tree, seqs)
            #print scores[node]
    #okay, now onto nodes that DO have children
    for node in tree:
        if node not in scores:
            scores = internalNode(node, scores, tree, subMatrix)

    return scores

def internalNode(node, scores, tree, subMatrix):
    score = {}
    #first, we want to make sure we score all the children
    for child in tree[node]:
        if child not in scores:
            scores = internalNode(child, scores, tree, subMatrix)
    choices = ['A', 'C', 'T', 'G'] #all the possible things value could be
    for choice in choices:
        low = 100000.0
        child1 = None
        #Let's us try to find the lowest score this can get from child1
        for entry in scores[tree[node][0]]:
            curr = (scores[tree[node][0]][entry][0] + subMatrix[choice][entry])
            if curr <= low:
                low = curr
                child1 = entry
        #print child1
        low2 = 100000.0
        child2 = None
        for entry in scores[tree[node][1]]:
            curr = (scores[tree[node][1]][entry][0] + subMatrix[choice][entry])
            if curr <= low2:
                low2 = curr
                child2 = entry
            #print child2

        score[choice] = (low+low2, child1, child2)
    scores[node] = score
    return scores

def leafNode(node, tree, seqs):
    #a function that just fills in the dict for a single leaf node
    score = {}
    value = seqs[node] #getting actual value at the node
    choices = ['A', 'C', 'T', 'G'] #all the possible things value could be
    for choice in choices:
        if choice == value:     #nice, we found the right one, this should cost 0
            score[choice] = (0.0, None, None)
        else: # fix, gives us the wrong output
            score[choice] = (100000.0, None, None)
    return score

def inferredStates(results,root):
    #Helper function to print the inferred sequences
    min = 100001 #arbitrary num over 100000
    result = ""
    minNode = ""
    rootnum = 0
    for choice in results:
        for node in results.get(choice):
            for num in results.get(choice).get(node):
                if num == 0:
                    minNode = node
                if num < min:
                    min = num
                    minNode = node
        if str(choice) == str(root):
            rootnum = min
        result += str(minNode)
        min = 100001
    return result,rootnum

 
                
            
def main():
     [tree, subMatrix, seqs, seqlen,root] = parseInput(argv)
     printStates(tree,subMatrix,seqs,seqlen,root)

def printStates(tree,subMatrix,seqs,seqlen,root):
    # Function prints out the inferred sequences
     print ""
     print "Inferred States:"
     tempdict = {}
     results = ""
     num = 0.0
     # loops through results
     for x in range(0,seqlen):
         for choice in seqs:
             nucleo = seqs.get(choice)
             tempdict[choice] = nucleo[x]
         seque, rootnum = inferredStates(wp(tree,subMatrix,tempdict),root)
         results += seque
         num += int(rootnum)
         #print rootnum
         tempdict = {}
     outcome = ""
     acc = 0
     
     for choices in tree:
         for y in range(0+acc,len(tree)*seqlen,len(tree)):
             outcome += results[y]
         print choices + ": " + outcome
         outcome = ""
         acc += 1

     print "score: " + str(num)
        

main()
