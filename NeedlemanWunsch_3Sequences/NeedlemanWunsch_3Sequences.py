#!usr/bin/python
import array
import numpy
  
# This function reads the two sequences from the file 'NeedlemanWunsch3_Input.txt'
# It also reads the standard distribution matrix from the file'pam250.txt'
def readInput():    
    # Read the sequences to be aligned
    Input_handler = open('NeedlemanWunsch3_Input.txt');
    lines = Input_handler.readlines();
    Input_handler.close();  
    
    # Process the sequences
    Sequence_a = lines[1].strip();
    Sequence_b = lines[3].strip();
    Sequence_c = lines[5].strip();
    print "Sequences to be aligned:"
    print Sequence_a
    print Sequence_b
    print Sequence_c
    print ""
    
    # Load a standard distribution matrix
    with open('pam250.txt') as f:
        Pam250 = [[x for x in ln.split()] for ln in f] 
        Pam250[0]=[0]+Pam250[0] 

    return Sequence_a, Sequence_b, Sequence_c, Pam250

# This function is used in the initialization of matrices
def g(k): return k*2

# This function is used to find the similarity score
def Sc(a, b, c, Pam250):
    return int(Pam250[Pam250[0].index(a)][Pam250[:][0].index(b)])+int(Pam250[Pam250[0].index(a)][Pam250[:][0].index(c)])+int(Pam250[Pam250[0].index(b)][Pam250[:][0].index(c)])

# This function is used to oversee the setup, population of the matrices required by the algorithm
def createMatrixes(Sequence_a, Sequence_b, Sequence_c, Pam250):

    # Create and initialise the 2D matrix for sequence 1 and sequence 2
    D_1_2 = twoDMatrix(Sequence_a, Sequence_b, Pam250);

    # Create and initialise the 2D matrix for sequence 2 and sequence 3
    D_2_3 = twoDMatrix(Sequence_b, Sequence_c, Pam250);

    # Create and initialise the 2D matrix for sequence 1 and sequence 3
    D_3_1 = twoDMatrix(Sequence_a, Sequence_c, Pam250);    
    
    # Create and initialise the 3D Alignment matrix D[i][j][k]
    D=numpy.zeros((len(Sequence_a)+1,len(Sequence_b)+1,len(Sequence_c)+1))
    
    size_i = len(D)
    size_j = len(D[0])
    size_k = len(D[0][0])
    for i in range(1, size_i):
        D[i][0][0]= g(i)
    
    for i in range(1, size_j):
        D[0][i][0]=g(i)
    
    for i in range(1, size_k):
        D[0][0][i]= g(i)

    for i in range(1, size_i):
        for j in range(1, size_j):
            Score = i + j
            D[i][j][0]=D_1_2[i][j]+g(Score)

    for j in range(1, size_j):
        for k in range(1, size_k):
            Score = j + k
            D[0][j][k]=D_2_3[j][k]+g(Score)
    
    for i in range(1, size_i):
        for k in range(1, size_k):
            Score = i + k
            D[i][0][k]=D_3_1[i][k]+g(Score)     

    # Populating the matrix
    for i in range(1, size_i):
        for j in range(1, size_j):
             for k in range(1, size_k):

                 Substitution = D[i-1][j-1][k-1]+Sc(Sequence_a[i-1],Sequence_b[j-1],Sequence_c[k-1], Pam250)
                 i_j_Gap = D[i-1][j-1][k]+Sc(Sequence_a[i-1],Sequence_b[j-1],'*', Pam250)
                 i_Gap_k = D[i-1][j][k-1]+Sc(Sequence_a[i-1],'*',Sequence_c[k-1], Pam250)
                 Gap_j_k = D[i][j-1][k-1]+Sc('*',Sequence_b[j-1],Sequence_c[k-1], Pam250)
                 i_Gap_Gap = D[i-1][j][k]+Sc(Sequence_a[i-1],'*','*', Pam250)
                 Gap_j_Gap = D[i][j-1][k]+Sc('*',Sequence_b[j-1],'*', Pam250)
                 Gap_Gap_k = D[i][j][k-1]+Sc('*','*',Sequence_c[k-1], Pam250)
                
                 D[i][j][k]= max(Substitution, i_j_Gap, i_Gap_k, Gap_j_k, i_Gap_Gap, Gap_j_Gap, Gap_Gap_k)
    return D, D[i][j][k] 

# This function has the logic to create two dimensional matrices
def twoDMatrix(Sequence_1, Sequence_2, Pam250):
    # Create the 2D Alignment matrix D[i][j]
    rows = len(Sequence_1) + 1;
    columns = len(Sequence_2) + 1;
    Gap = -8;
    D = [[0 for x in range(columns)] for x in range(rows)]
    
    # Initialize the matrix
    D[0][0] = 0 
    for i in range(1,rows):
        D[i][0] = D[i-1][0] + int(Gap)
    for j in range(1,columns):
        D[0][j] = D[0][j-1] + int(Gap)
    
    # Populate the matrix
    for i in range(1, rows):
        for j in range(1, columns):
            Substitution= D[i-1][j-1]+int(Pam250[Pam250[0].index(Sequence_1[i-1])][Pam250[:][0].index(Sequence_2[j-1])]) 
            Insertion = D[i][j-1] + Gap
            Deletion = D[i-1][j] + Gap
            D[i][j] = max([Substitution, Insertion, Deletion])        
    return D

# This function performs the traceback operation 
def traceBack(D, Sequence_a, Sequence_b, Sequence_c, Pam250):

    Statement_1 = ""
    Statement_2 = ""       
    Statement_3 = ""
    
    i = len(Sequence_a) + 1
    j = len(Sequence_b) + 1
    k = len(Sequence_c) + 1

    i = i-1
    j = j-1
    k = k-1
    
    while True:
        
        if i==0 and j==0 and k ==0:
            break;
        if i==0 or j==0 or k==0:
            # When i is 0
            if i==0 and j!=0 and k!=0:
                Statement_1 = '-' + Statement_1
                Statement_2 = Sequence_b[j] + Statement_2
                Statement_3 = Sequence_c[k] + Statement_3
                j = j-1
                k = k-1
            # When j is 0
            elif i!=0 and j==0 and k!=0:
                Statement_1 = Sequence_a[i] + Statement_1
                Statement_2 = '-' + Statement_2
                Statement_3 = Sequence_c[k] + Statement_3
                i = i-1
                k = k-1 
            # When k is 0
            elif i!=0 and j!=0 and k==0:
                Statement_1 = Sequence_a[i] + Statement_1
                Statement_2 = Sequence_b[j] + Statement_2
                Statement_3 = '-' + Statement_3
                i = i-1
                j = j-1
            # When i and j are 0
            elif i==0 and j==0 and k!=0:
                Statement_1 = '-' + Statement_1
                Statement_2 = '-' + Statement_2
                Statement_3 = Sequence_c[k] + Statement_3
                k = k-1
            # When i and k are 0
            elif i==0 and j!=0 and k==0:
                Statement_1 = '-' + Statement_2
                Statement_2 = Sequence_c[j] + Statement_2
                Statement_3 = '-' + Statement_3
                j = j-1  
            # When j and k are 0
            elif i!=0 and j==0 and k==0:
                Statement_1 = Sequence_a[i] + Statement_1
                Statement_2 = '-' + Statement_2
                Statement_3 = '-' + Statement_3
                i = i-1
        else:
            if D[i][j][k] == D[i-1][j-1][k-1] + Sc(Sequence_a[i-1], Sequence_b[j-1], Sequence_c[k-1], Pam250):
                Statement_1 = Sequence_a[i-1] + Statement_1
                Statement_2 = Sequence_b[j-1] + Statement_2
                Statement_3 = Sequence_c[k-1] + Statement_3
                i = i-1
                j = j-1
                k = k-1
            elif D[i][j][k] == D[i-1][j-1][k] + Sc(Sequence_a[i-1], Sequence_b[j-1], '-', Pam250):
                Statement_1 = Sequence_a[i-1] + Statement_1
                Statement_2 = Sequence_b[j-1] + Statement_2
                Statement_3 = '-' + Statement_3
                i = i-1
                j = j-1
            elif D[i][j][k] == D[i-1][j][k-1] + Sc(Sequence_a[i-1], '-', Sequence_c[k-1], Pam250):
                Statement_1 = Sequence_a[i-1] + Statement_1
                Statement_2 = '-' + Statement_2
                Statement_3 = Sequence_b[k-1] + Statement_3
                i = i-1
                k = k-1
            elif D[i][j][k] == D[i][j-1][k-1] + Sc('-', Sequence_b[j-1], Sequence_c[k-1], Pam250):
                Statement_1 = '-' + Statement_1
                Statement_2 = Sequence_b[j-1] + Statement_2
                Statement_3 = Sequence_c[k-1] + Statement_3
                j = j-1
                k = k-1
            elif D[i][j][k] == D[i-1][j][k] + Sc(Sequence_a[i-1], '-', '-', Pam250):
                Statement_1 = Sequence_a[i-1] + Statement_1
                Statement_2 = '-' + Statement_2
                Statement_3 = '-' + Statement_3
                i = i-1
            elif D[i][j][k] == D[i][j-1][k] + Sc('-', Sequence_b[j-1], '-', Pam250):
                Statement_1 = '-' + Statement_1
                Statement_2 = Sequence_b[j-1] + Statement_2
                Statement_3 = '-' + Statement_3
                j = j-1
            elif D[i][j][k] == D[i][j][k-1] + Sc('-', '-', Sequence_c[k-1], Pam250):
                Statement_1 = '-' + Statement_1
                Statement_2 = '-' + Statement_2
                Statement_3 = Sequence_c[k-1] + Statement_3
                k = k-1
    return Statement_1, Statement_2, Statement_3
                
# alignSequences does the task of three 2D matrices and one 3D matrix
# It also performs the traceback operation
def alignSequences():
    
    # Read the sequences
    Sequence_a, Sequence_b, Sequence_c, Pam250 = readInput()
    rows = len(Sequence_a) + 1;
    columns = len(Sequence_b) + 1;
    height = len(Sequence_c) + 1;
    if (str.isalpha(Sequence_a) and str.isalpha(Sequence_b) and str.isalpha(Sequence_c)):              
        D, Score = createMatrixes(Sequence_a, Sequence_b, Sequence_c, Pam250);  
        Statement_1, Statement_2, Statement_3 = traceBack(D, Sequence_a, Sequence_b, Sequence_c, Pam250);
        print "*******************************************************************************" 
        print "Score of the alignment:",int(Score)
        print ""
        print 'Optimal multiple sequence alignment:'
        print Statement_1
        print Statement_2
        print Statement_3
        print "*******************************************************************************" 
    else:
         print "The input sequences are invalid"
def main():
    alignSequences();

if __name__ == "__main__":
    main()  

   

