#!usr/bin/python
import array
import numpy

# In the WPGMA method, we can align multiple sequences.
# As input, we require more than 2 sequences and a distance matrix

# This function reads the distance matrix from the file 'WPGMA_Input.txt'
def readInput():      
    matrix = []
    with open('WPGMA_Input.txt') as f:
        matrix = [[x for x in ln.split()] for ln in f]       
    matrix = numpy.asarray(matrix)
    matrix = matrix.astype(numpy.float)
    length = len(matrix)
    print "Distance matrix:"
    print matrix
    print ""   
    print "Number of sequences:", length
    print ""
    return matrix, length

# This function finds the minimum value in the distance matrix
def matrixMinimum(matrix, length):
    min_index_i = 0
    min_index_j = 0
    minimum=float('inf')
    for i in xrange (length):

        temp = matrix[i]
        min_value = numpy.min(temp[numpy.nonzero(temp)])
        j = temp.tolist().index(min_value)

        if min_value < minimum: 
            min_index_i = i
            min_index_j = j
            minimum = min_value
   
    return min_index_i,min_index_j

def wpgma(matrix, length,dictionary):
    leaves = []
    count=0
    for i in xrange (0, length):    
        leaves.append("S"+str(i+1))
    
    numberOfClusters=i+1
    
    while(length>1):
        
        
        numberOfClusters=numberOfClusters+1
        count=count+1
        min_index_i,min_index_j = matrixMinimum(matrix, length)
        
        leaves.append("S"+str(numberOfClusters))  
        distance=matrix[min_index_i][min_index_j]/float(2)
        
        size=0
        if leaves[min_index_i] not in dictionary.keys():
            size1=1
            distance1=distance
        else:
            size1=dictionary[leaves[min_index_i]][4]
            distance1=distance-max(dictionary[leaves[min_index_i]][0],dictionary[leaves[min_index_i]][2])
        
        if leaves[min_index_j] not in dictionary.keys():
            size2=1
            distance2=distance
        else:
            size2=dictionary[leaves[min_index_j]][4]
            distance2=distance-max(dictionary[leaves[min_index_j]][0],dictionary[leaves[min_index_j]][2])
        dictionary["S"+str(numberOfClusters)]=[distance1,leaves[min_index_i],distance2,leaves[min_index_j],size1+size2]
        
        # Create a new row and column
        matrix = numpy.insert(matrix, length, values=float(0), axis=0)
        matrix = numpy.insert(matrix, length, values=float(0), axis=1)
        
        for i in xrange (0, length):
            matrix[-1][i]=matrix[i][-1] = (size1*matrix[i][min_index_i] + size2*matrix[i][min_index_j])/float(size1+size2)
        
        # Delete the minimum value
        if min_index_i < min_index_j:
            matrix = numpy.delete(matrix, min_index_i, 0)
            matrix = numpy.delete(matrix, min_index_i, 1)
            matrix = numpy.delete(matrix, (min_index_j)-1, 0)
            matrix = numpy.delete(matrix, (min_index_j)-1, 1)
            length = len(matrix)
            del leaves[min_index_j]
            del leaves[min_index_i]

        else:
            matrix = numpy.delete(matrix, min_index_i, 0)
            matrix = numpy.delete(matrix, min_index_i, 1)
            matrix = numpy.delete(matrix, min_index_j, 0)
            matrix = numpy.delete(matrix, min_index_j, 1)            
            length = len(matrix)
            del leaves[min_index_i]
            del leaves[min_index_j]
        
    return "S"+str(numberOfClusters)

def printCluster(dictionary,finalCluster):
    stack=[]
    result=[]
    stack.append(finalCluster)
    while stack:
        
        current=stack.pop()
        if isinstance( current , float ):
            if isinstance( current_prev , float ):
                result.pop()
                result.append(")")
            result.append(":"+str(current))
            result.append(",")
            
        elif current in dictionary.keys():
        
            stack.append(dictionary[current][0])
            stack.append(dictionary[current][1])
            stack.append(dictionary[current][2])
            stack.append(dictionary[current][3])
            result.append("(")
        else:
            result.append(current)
        current_prev=current
        
    result.pop()
    result.append(")")  
    return result

if __name__ == "__main__":
    matrix, length = readInput()
    dictionary={}
    finalCluster = wpgma(matrix, length,dictionary)
    result=printCluster(dictionary,finalCluster)
    result=''.join(result)
 
    result=result+";"
    from ete2 import Tree
    tree=Tree(result)
    print "*******************************************************************************" 
    print "WPGMA Resultant Clustering:"
    print ""   
    print result 
    print ""    
    print tree
    print "*******************************************************************************" 
