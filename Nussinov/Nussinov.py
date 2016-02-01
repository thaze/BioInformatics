#!usr/bin/python
import array
import numpy

'''
This function reads the aligned sequences from the file 'Nussinov_Input.txt'
'''
def readInput():
    # Read the sequences to be aligned
    Input_handler = open('Nussinov_Input.txt');
    lines = Input_handler.readlines();
    Input_handler.close();  
    
    # Process the sequences
    Sequence = lines[0].strip();
    print "Input Sequence:"
    print Sequence
    print ""

    # Calculate the length of the sequences
    length = len(Sequence);
    return Sequence, length;

'''
This function creates a matrix for the sequence.
'''
def createMatrix(Sequence, length):
    matrix = [[0 for x in range (length)] for x in range (length)]
    matrix = numpy.asarray(matrix)
    i=0
    j=1
    k=1
    comparison = 0
    
    for index in range (length):
        i=0
        j=k
        while j<length:
            comparison = 0
            if (Sequence[i] == 'C' and Sequence[j] == 'G') or (Sequence[i] == 'G' and Sequence[j] == 'C'):
                comparison = 1
            elif (Sequence[i] == 'A' and Sequence[j] == 'U') or (Sequence[i] == 'U' and Sequence[j] == 'A'):
                comparison = 1
            max_number=[]
            max_value=0
            for scope in range(i+1,j):
                    max_number.append(matrix[i][scope]+matrix[scope+1][j])
                    max_value=max(max_number)
            matrix[i][j]=max(matrix[i+1][j],matrix[i][j-1],matrix[i+1][j-1]+comparison,max_value)
            
            i=i+1
            j=j+1
        k=k+1
    print "Maximal number of base pairs", matrix[0][len(Sequence)-1]
    print ""
    return matrix

'''
This function performs the traceback operation to find the base pairs
'''
def traceBack(matrix, Sequence, length):
    # Initialisation
    coordinates=[]
    result=["."]*length
    stack = []
    comparison = 0
    stack.append([0, length-1])
        
    while stack:
        i=stack[-1][0]
        j=stack[-1][1]
        if Sequence[i]=='C' and Sequence[j]=='G':
            comparison=1        
        if Sequence[i]=='G' and Sequence[j]=='C':
            comparison=1        
        if Sequence[i]=='A' and Sequence[j]=='U':
            comparison=1
        if Sequence[i]=='U' and Sequence[j]=='A':
            comparison=1

        stack.pop()
        if i>=j:
           continue
        elif matrix[i+1,j]==matrix[i,j]:
            stack.append([i+1,j])
        elif matrix[i,j-1]==matrix[i,j]:
            stack.append([i,j-1])
        elif matrix[i+1,j-1]+comparison==matrix[i,j]:
            coordinates.append([i,j])
            result[i]="("
            result[j]=")"  
            stack.append([i+1,j-1])
        else:
            for k in range(i+1,j):
                if matrix[i,k]+matrix[k+1,j]==matrix[i,j]:
                    stack.append([k+1,j])
                    stack.append([i,k])
                    break
    print 'Base Pair representation:'
    print Sequence
    print ""
    print "Optimal structure dot-bracket:"
    print ''.join(result)
    print ""
    print "Optimal structure coordinates:"
    print coordinates
    return result    
    
'''
This function invokes the functions 'readInput' , 'createMatrix' and 'traceBack'
'''
def findBasePairs():
    Sequence, length = readInput();
    if (str.isalpha(Sequence)):
        matrix = createMatrix(Sequence, length);
        result = traceBack(matrix, Sequence, length);        
    else:
         print "Error: The input sequence is invalid"

'''
The main function
'''
def main():
    findBasePairs();

if __name__ == "__main__":
    main()  
