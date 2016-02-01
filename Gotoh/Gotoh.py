#!usr/bin/python
import array
import numpy


'''
This class defines the nodes in the linked list. Each node has 3 data values and a pointer to the next node. 
'''
class Node(object):
 
    def __init__(self, data1, data2, data3, next):
        self.data1 = data1
        self.data2 = data2
        self.data3 = data3
        self.next = next

'''
A linked list is used during traceback.
The LinkedList class is used to add and show nodes, delete nodes till a mentioned value and print the alignment
'''
class LinkedList(object):
    head = None
    tail = None

    # This function can be used to display all the nodes in the linked list.
    def show(self):
        current_node = self.head
        while current_node is not None:
            current_node = current_node.next
        print None

    # This function is invoked with the three data items required to create a new node in the linked list.           
    def add(self, data1, data2, data3):
        node = Node(data1, data2, data3, None)
        if self.head is None:
            self.head = self.tail = node
        else:
            self.tail.next = node
        self.tail = node
               
    # This function aids in the process of acquiring multiple optimal alignments. Following the print of an optimal 
    # alignment, this function is called to delete the nodes corresponding to the section of the path already printed.
    # If traceback started at D[3][5] and there was a split at D[2][4], then the nodes prior to it will be deleted.
    def delete_till(self, node_value1, node_value2, node_value3):
        current_node = self.head
        previous_node = None

        while current_node is not None:
            if current_node.data1 == node_value1 and current_node.data2 == node_value2 and current_node.data3 == node_value3:
                current_node.next=None
                self.tail = current_node
                
            # The next iteration
            previous_node = current_node
            current_node = current_node.next
    
    # This function is invoked when the alignment needs to be displayed               
    def displayList(self, Sequence_a, Sequence_b):  
        Statement_1 = ''
        Statement_2 = ''
        current_head = self.head
        temp1 = self.head
        temp2 = self.head
        temp2 = temp2.next
        
        while temp2:
            
            if (temp2.data1 == (temp1.data1)-1) and (temp2.data2 == (temp1.data2)-1):
                Statement_1 =  Sequence_a[-1] + Statement_1
                Statement_2 =  Sequence_b[-1] + Statement_2
                Sequence_a=Sequence_a[:-1]
                Sequence_b=Sequence_b[:-1]
            elif (temp2.data1 == ((temp1.data1)-1)) and (temp2.data2 == temp1.data2):
                Statement_1 = Sequence_a[-1] + Statement_1
                Statement_2 = '-' + Statement_2
                Sequence_a=Sequence_a[:-1]
            elif (temp2.data1 == temp1.data1) and (temp2.data2 == (temp1.data2)-1):
                Statement_1 = '-' + Statement_1
                Statement_2 = Sequence_b[-1] + Statement_2
                Sequence_b=Sequence_b[:-1]
                
            temp1 = temp1.next
            temp2 = temp2.next
            
 
            
        return Statement_1, Statement_2 

'''
This function reads the two sequences, alpha and beta from the file 'Gotoh_Input.txt'
It also reads the standard distribution matrix from the file'pam250.txt'
'''
def readSequences():    
    # Read the sequences from a file
    file_handler = open('Gotoh_Input.txt')
    lines = file_handler.readlines()
    file_handler.close()
    
    # Process the sequences
    Sequence_a = lines[1].strip();
    Sequence_b = lines[3].strip();
    alpha = lines[5].strip();
    alpha = int(alpha)
    beta = lines[7].strip();
    beta = int(beta)
    print "alpha:",alpha  
    print "beta:",beta
    print ""
    print "Sequences to be aligned:"
    print Sequence_a
    print Sequence_b
    print ""

    # Calculate the length of the sequences
    len_a = len(Sequence_a)
    len_b = len(Sequence_b)

    # Load a standard distribution matrix
    with open('pam250.txt') as f:                         
     Pam250 = [[x for x in ln.split()] for ln in f]
     Pam250[0]=[0]+Pam250[0] 
      
    return Sequence_a, Sequence_b, alpha, beta, len_a, len_b, Pam250

'''
alignSequences does the task of creating, initializing and populating the D[i][j], P[i][j] and Q[i][j] matrices
It also performs the traceback operation
'''
def AlignSequences():
    # Initialize the symbol string 
    Symbol_String = ""
    Statement_1=""
    Statement_2=""    
    myList = LinkedList()
    
    # Store_Cells is used to remember the cells visited and the one's where the path splits up
    Store_Cells = [] 
    
    Sequence_a, Sequence_b, alpha, beta, len_a, len_b, Pam250 = readSequences()
    Infinity = float("inf")

    # Calculate the number rows and columns in the matrix
    rows = len_a + 1;
    columns= len_b + 1;

    if (str.isalpha(Sequence_a) and str.isalpha(Sequence_b)):
        # Create and Initialize the Pi,j matrix
        P = [[0 for x in range(columns)] for x in range(rows)]
        P[0][0] = 0
        for i in range(1,rows):
            P[i][0] = '-'; 
        for j in range(1,columns):
            P[0][j] = -Infinity;
            
        # Create and Initialize the Qi,j matrix
        Q = [[0 for x in range(columns)] for x in range(rows)]
        Q[0][0] = 0
        for i in range(1,rows):
            Q[i][0] = -Infinity; 
        for j in range(1,columns):
            Q[0][j] = '-';
    
        # Create and Initialize the Di,j matrix
        D = [[0 for x in range(columns)] for x in range(rows)]
        D[0][0] = 0
        for i in range(1,rows):
            D[i][0] = alpha + (beta * i); 
        
        for j in range(1,columns):
            D[0][j] = alpha + (beta * j); 
            
        # Populate the matrices         
        for i in range(1, rows):
            for j in range(1, columns):
               
                # P Matrix
                case1 = D[i-1][j] + (alpha + (beta * 1))
                case2 = P[i-1][j] + beta
                P[i][j] = max([case1, case2])
                
                # Q Matrix
                case3 = D[i][j-1] + (alpha + (beta * 1))
                case4 = Q[i][j-1] + beta
                Q[i][j] = max([case3, case4])
    
                # D matrix
                D[i][j] = D[i-1][j-1]+int(Pam250[Pam250[0].index(Sequence_a[i-1])][Pam250[:][0].index(Sequence_b[j-1])])
                D[i][j] = max([D[i][j], P[i][j], Q[i][j]])
        print "Score:",D[i][j]                      
            
        #print "Now that we have the aligment matrix populated with the scores, lets proceed with traceback"    
        #Traceback starts at the cell with the highest value
        i = len_a   
        j = len_b
    
        Store_Cells.append([0,0,i,j,'D'])
       
        # This while loop holds good as long as Store_Cells is not empty
        while Store_Cells: 
            i = Store_Cells[-1][2]   
            j = Store_Cells[-1][3]
            M = Store_Cells[-1][4]
            myList.add(i,j,M)
            
            Store_Cells.pop()
    
            if i==0 and j==0:

                Statement_1, Statement_2 = myList.displayList(Sequence_a, Sequence_b)
                print "*******************************************************************************" 
                print "Pairwise Alignment(s):" 
                print ""
                print Statement_1
                print Statement_2
                print "*******************************************************************************" 
                
                if Store_Cells:
                # After printing the alignment, we can delete the nodes from the list until p and q are reached
                    p = Store_Cells[-1][0]
                    q = Store_Cells[-1][1]
                    M = Store_Cells[-1][4]
                    
                    myList.delete_till(p,q,M)
                continue
                
            if i==0 or j==0 and M=='D':
                if i==0:
                    Store_Cells.append([0,j,0,j-1,'D'])
                else:
                    Store_Cells.append([i,0,i-1,0,'D'])
                continue
        
            # Check whether i and j resulted from substitution, insertion or deletion.
            if M == 'D':
                if D[i][j] == D[i-1][j-1]+int(Pam250[Pam250[0].index(Sequence_a[i-1])][Pam250[:][0].index(Sequence_b[j-1])]):
                    Store_Cells.append([i,j,i-1,j-1, 'D'])
    
                if D[i][j] == P[i][j]:
                    Store_Cells.append([i,j,i,j, 'P'])
                
                if D[i][j] == Q[i][j]:
                    Store_Cells.append([i,j,i,j,'Q'])
                    
            elif M=='P':
                # If the value resulted from deletion
                if P[i][j] == D[i-1][j] + (alpha + (beta * 1)):
                    Store_Cells.append([i,j,i-1,j,'D'])
                
                # If the value resulted from insertion
                if P[i][j] == P[i-1][j] + beta:
                    Store_Cells.append([i,j,i-1,j, 'P'])
            else:
                # Q Matrix
                if Q[i][j] == D[i][j-1] + (alpha + (beta * 1)):
                    Store_Cells.append([i,j,i,j-1, 'D'])
                
                if Q[i][j]== Q[i][j-1] + beta:
                    Store_Cells.append([i,j,i,j-1, 'D'])
    else:
         print "Error: The input sequences are invalid"

def main():
    AlignSequences();
    
if __name__ == "__main__":
    main()    
        

