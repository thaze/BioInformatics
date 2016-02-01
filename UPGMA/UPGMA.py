#!usr/bin/python
import array
import numpy

# In the UPGMA method, we can align multiple sequences.
# As input, we require more than 2 sequences and a distance matrix

# This function reads the distance matrix from the file 'UPGMA_Input.txt'
def readInput():      
    matrix = []
    with open('UPGMA_Input.txt') as f:
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

# This function locates the smallest entry in the distance matrix.
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

'''
This function deals with the processing involved in the clustering according to the upgma algorithm. 
The dictionary data structure is used.
Initially, leaves (tree nodes) need to be created for each sequence.
These leaves are stored in an array under the same name 'leaves'.
The main loop works under the condition that there is more than one sequence. 
Within the loop, the smallest value in the input matrix is retrieved. 
If this smallest value i.e. min_index_i is not a key value in the dictionary, we set the size as '1' and the distance is calculated 
using the upgma formula. On the otherhand, if it is a key value, 

If min_index_j is not a key value in the dictionary


Next, a new row and column in the matrix is created for the new leaf and the distance values are calculated
The row and column with the minimum distance value is deleted.


'''

def upgma(matrix, length,dictionary):
    leaves = []
    count = 0
    
    for i in xrange (0, length):    
        leaves.append("S"+str(i+1))
    
    numberOfClusters = i+1
    
    while(length>1):
        
        numberOfClusters = numberOfClusters+1
        count = count+1

        min_index_i,min_index_j = matrixMinimum(matrix, length)
        
        leaves.append("S"+str(numberOfClusters))  
        distance = matrix[min_index_i][min_index_j]/float(2)
        
        size = 0
        if leaves[min_index_i] not in dictionary.keys():
            size = 1
            distance1 = distance
        else:
            size = dictionary[leaves[min_index_i]][4]
            distance1 = distance-max(dictionary[leaves[min_index_i]][0],dictionary[leaves[min_index_i]][2])
        
        if leaves[min_index_j] not in dictionary.keys():
            size = size+1
            distance2 = distance
        else:
            size = size+dictionary[leaves[min_index_j]][4]
            distance2 = distance-max(dictionary[leaves[min_index_j]][0],dictionary[leaves[min_index_j]][2])
        dictionary["S"+str(numberOfClusters)] = [distance1,leaves[min_index_i],distance2,leaves[min_index_j],size]
        
        # Create a new row and column
        matrix = numpy.insert(matrix, length, values=float(0), axis=0)
        matrix = numpy.insert(matrix, length, values=float(0), axis=1)

        for i in xrange (0, length):
            matrix[-1][i]=matrix[i][-1] = (matrix[i][min_index_i] + matrix[i][min_index_j])/2
        
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

'''
This is the function used to print the cluster resulting as a result of the upgma method
'''

def printCluster(dictionary,finalCluster):
    stack = []
    result = []
    stack.append(finalCluster)
    while stack:
        
        current = stack.pop()
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
        current_prev = current
        
    result.pop()
    result.append(")")  
    return result

if __name__ == "__main__":
    matrix, length = readInput()
    dictionary={}
    finalCluster = upgma(matrix, length,dictionary)
    result=printCluster(dictionary,finalCluster)
    result=''.join(result)
    result=result+";"
    
    from ete2 import Tree
    tree=Tree(result)
    print "*******************************************************************************" 
    print "UPGMA Resultant Clustering:"
    print ""   
    print result 
    print ""    
    print tree
    print "*******************************************************************************" 
