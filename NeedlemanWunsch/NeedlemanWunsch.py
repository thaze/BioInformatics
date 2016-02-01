#!usr/bin/python
import array
import numpy

'''
This class defines the nodes in the linked list. Each node has two data values and a pointer to the next node. 
'''
class Node(object):
 
    def __init__(self, data1, data2, next):
        self.data1 = data1
        self.data2 = data2
        self.next = next

'''
A linked list is used during traceback.
The LinkedList class handles the functionality to add and show nodes, delete nodes until a specific value and print 
the entire alignment. 
'''
class LinkedList(object):
    head = None
    tail = None

    # This function can be used to display all the nodes in the linked list.
    def show(self):
        print "The nodes are:"
        current_node = self.head
        while current_node is not None:
            print current_node.data1, current_node.data2
            current_node = current_node.next
        print None

    # This function is invoked with the two data items required to create a new node in the linked list.    
    def add(self, data1, data2):
        node = Node(data1, data2, None)
        if self.head is None:
            self.head = self.tail = node
        else:
            self.tail.next = node
        self.tail = node
       
    # This function aids in the process of acquiring multiple optimal alignments. Following the print of an optimal 
    # alignment, this function is called to delete the nodes corresponding to the section of the path already printed.
    # If traceback started at D[3][5] and there was a split at D[2][4], then the nodes prior to it will be deleted.
    def delete_till(self, node_value1, node_value2):
        current_node = self.head
        previous_node = None
        while current_node is not None:
            if current_node.data1 == node_value1 and current_node.data2 == node_value2:
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
This function reads the two sequences from the file 'NeedlemanWunsch_Input.txt'
It also reads the standard distribution matrix from the file'PAM250.txt'
'''
def readInput():    
    # Read the sequences to be aligned
    Input_handler = open('NeedlemanWunsch_Input.txt');
    lines = Input_handler.readlines();
    Input_handler.close();  
    
    # Process the sequences
    Sequence_a = lines[1].strip();
    Sequence_b = lines[3].strip();   
    print "Sequences to be aligned:"
    print Sequence_a
    print Sequence_b
    print ""
    
    # Calculate the length of the sequences
    len_a = len(Sequence_a);
    len_b = len(Sequence_b); 
    Gap = -8;    
    print "Gap:", Gap
    print ""

    # Load a standard distribution matrix
    with open('pam250.txt') as f:                         
     Pam250 = [[x for x in ln.split()] for ln in f]
     Pam250[0]=[0]+Pam250[0]

    return Sequence_a, Sequence_b, Gap, len_a, len_b, Pam250

'''
alignSequences does the task of creating, initializing and populating the D[i][j] matrix
It also performs the traceback operation
'''
def alignSequences():
    # Initialize the symbol string 
    Symbol_String = ""
    Statement_1=""
    Statement_2=""    
    myList = LinkedList()
    
    # Store_Cells is used to remember the cells visited and the one's where the path splits up
    Store_Cells = []    
    # Retrieve the sequences    
    Sequence_a, Sequence_b, Gap, len_a, len_b, Pam250 = readInput()

    # Calculate the number rows and columns in the matrix based on the length of the sequencs
    rows = len_a + 1;
    columns= len_b + 1;
    if (str.isalpha(Sequence_a) and str.isalpha(Sequence_b)):    
        # Create the 2D Alignment matrix D[i][j]
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
                Substitution= D[i-1][j-1]+int(Pam250[Pam250[0].index(Sequence_a[i-1])][Pam250[:][0].index(Sequence_b[j-1])]) 
                Insertion = D[i][j-1] + Gap
                Deletion = D[i-1][j] + Gap
                D[i][j] = max([Substitution, Insertion, Deletion])        
        print "The score of the alignment is:", D[i][j]
          
        #Traceback starts at the cell D[i,j]
        Store_Cells.append([0,0,i,j])
        print "*******************************************************************************" 
        print "Pairwise Alignment(s):" 
        
        # This while loop holds good as long as Store_Cells is not empty
        while Store_Cells:
        
            i = Store_Cells[-1][2]   
            j = Store_Cells[-1][3]
    
            myList.add(i,j)
            Store_Cells.pop()
    
            if i==0 and j==0:
                Statement_1, Statement_2 = myList.displayList(Sequence_a, Sequence_b)                     
                print Statement_1
                print Statement_2
                print ""        
                
                if Store_Cells:
                # After printing the alignment, we can delete the nodes from the list until p and q are reached
                    p = Store_Cells[-1][0]
                    q = Store_Cells[-1][1]
                    myList.delete_till(p,q)
                continue
            
            if i==0 or j==0:
                if i==0:
                    Store_Cells.append([0,j,0,j-1])
                else:
                    Store_Cells.append([i,0,i-1,0])
                continue
                  
            
            # Check whether i and j resulted from substitution, insertion or deletion.
            # If the value resulted from substitution 
            if D[i][j] == D[i-1][j-1]+int(Pam250[Pam250[0].index(Sequence_a[i-1])][Pam250[:][0].index(Sequence_b[j-1])]):
                Store_Cells.append([i,j,i-1,j-1])
            
            # If the value resulted from deletion
            if D[i][j] == D[i-1][j] + Gap:
                Store_Cells.append([i,j,i-1,j])
                
            # If the value resulted from insertion
            if D[i][j] == D[i][j-1] + Gap:
                Store_Cells.append([i,j,i,j-1])
        print "*******************************************************************************" 
    else:
         print "Error: The input sequences are invalid"                
    
def main():
    alignSequences();

if __name__ == "__main__":
    main()  

   

