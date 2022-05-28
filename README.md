# Code source Structure 
The Program consists of 3 parts: the Main **(Hirschberg.py)**, and 2 modules: the hb  **(Hirschberg) class (hb.py)**, and the **nmw (Needleman-Wunsch) class (nmw.py)**

- **The hb class:** For the creation of “hb” objects, after which we call its Hirschberg method  to create the list of aligned sequences. 
- **The nmw class:** For the creation of “nmw” objects, after which we call its  EnumerateAlignments method to create the list of aligned sequences.
- **The main:** The starting point of the program’s execution. Results in both the Needleman Wunsch and Hirschberg Alignment.

## The Main (Hirschberg.py)
**The program takes 8 arguments:**

 1. 03 optional arguments: 
- “-t” to Print i & j iterations in the Hirschberg function. 
- “-f” to indicate if the input sequences are files. 
- “-l” to indicate if the file input sequences are aligned by lines & not characters. 
 2. 05 positional arguments: 
- “gap”: an integer parameter = the gap penalty “g” for the construction of F Score Matrix.
- “match”: an integer parameter = the match score “m” for the construction of F Score  Matrix". 
- “diff”: an integer parameter = the diff score “d” for the construction of F Score Matrix. ▪ “A”: The 1st sequence or file to be aligned. 
- “B”: The 2nd sequence or file to be aligned
- We create the ArgumentParser object to create the arguments and then get their values.
```python
import argparse 
import nmw 
import hb 
parser = argparse.ArgumentParser() 
# Creating optional arguments 
parser.add_argument("-t", action="store_true", help="Print i & j iterations") parser.add_argument("-f", action="store_true", help="The input is a file") parser.add_argument("-l", action="store_true", help="The file input sequences  are lines & not characters") 
# Creating positional arguments 
parser.add_argument("gap", type=int, help="Parameter gap penalty \"g\" for  the construction of F Score Matrix") 
parser.add_argument("match", type=int, help="Parameter match \"m\" for the  construction of F Score Matrix") 
parser.add_argument("diff", type=int, help="Parameter diff \"d\" for the  construction of F Score Matrix") 
parser.add_argument("A", help="The 1st sequence or file to be aligned") parser.add_argument("B", help="The 2nd sequence or file to be aligned") # Getting Arguments 
args = parser.parse_args() 
g = args.gap 
m = args.match 
d = args.diff 
t = args.t 
l = args.l 
f = args.f
```
-  Then the “hb” and “nmw” objects are instantiated. 
```python
align_nmw = nmw.nmw(g, m, d) 
align_hb = hb.hb(g, m, d)
```
- ‘f’ = True: 
    1) ‘l’ = True: 
    ```python
    if args.f: 
      with open(args.A, 'r') as A, open(args.B, 'r') as B: 
        if l: # a and b are lists of lines 
          a = A.read().splitlines() 
          b = B.read().splitlines() 
          print("\nAlignments using the NEEDLEMAN-WUNSCH Alignment ----\n")  f = align_nmw.F(a, b) 
          w = [] 
          z = [] 
          ww, zz = align_nmw.EnumerateAlignments_lines(a, b, f, w, z)  for w, z in zip(ww, zz): 
          for i, j in zip(w, z): 
            if i == j: 
              print("=", i, "\n=", j) 
        else: 
          print("<", i, "\n>", j) 
          print("\nAlignments using the HIRSCHBERG Alignment ----------\n")  ww, zz = align_hb.Hirschberg_lines(a, b) 
          for w, z in zip(ww, zz): 
          for i, j in zip(w, z): 
          if i == j: 
            print("=", i, "\n=", j) 
          else: 
            print("<", i, "\n>", j)
    ```
    2) ‘l’ = False: 
    ```python
     else: # a and b are strings 
        a = A.read().splitlines() 
        b = B.read().splitlines() 
        a = ''.join(a) 
        b = ''.join(b) 
        print("\nAlignments using the NEEDLEMAN-WUNSCH Alignment ----\n")  f = align_nmw.F(a, b) 
        w = "" 
        z = "" 
        ww, zz = align_nmw.EnumerateAlignments(a, b, f, w, z)  
        for w, z in zip(ww, zz):
            print(w) 
            print(z) 
        print("\nAlignments using the HIRSCHBERG Alignment ----------\n")  ww, zz = align_hb.Hirschberg(a, b) 
        alignments = zip(ww, zz) 
        alignments = list(set(alignments)) # Eliminates duplicate  Alignments 
        for w, z in alignments: 
            print(w) 
            print(z) 
    ```
**N.B** We use the list(set(alignments)) to eliminate the possible duplicate results of the Hirschberg method.
    
## Needleman-Wunsch Class (nmw.py)
When an object is instantiated from the nmw class, the “g”, “d”, and “m” parameters for the creation of the “f” score matrix are initialized, as well as the lists “ww” and “zz” where the aligned sequences of “a” and “b” respectively will be stored.
```python
class nmw: 

def __init__(self, g, m, d):
    self.g = g
    self.d = d
    self.m = m
    self.ww = []
    self.zz = []
```
**Methods:**
- Compare (x, y): used in the creation of “f” and in the back-tracing of the “f” elements origins, it compares x and y and returns m (reward) if match and d (penalty) if mismatch.
```python
def compare(self, x, y):  # Compares sequence elements x and y to determine the match/difference score
        if x == y:
            return self.m
        else:
            return self.d
```
-	F (a, b): creates the score matrix “f” of sequences “a” and “b”. The matrix is created by initializing the first line and column by multiplying their indices in the gap penalty (as moving right or down would result in a gap), then the rest of the matrix is filled using the Bellman's principle of optimality (calculating the step that results in the best score; moving diagonally means a match “m” or a mismatch “d”, moving right or down means a gap “g”).
```python
def F(self, a, b):  # Creates & returns the score matrix f
        r = len(a)  # number of rows
        c = len(b)  # number of columns
        f = [[0 for _ in range(c + 1)] for _ in range(r + 1)]  # a matrix of zeros
        for i in range(r + 1):  # filling the first column
            f[i][0] = i * self.g
        for j in range(c + 1):  # filling the first row
            f[0][j] = j * self.g
        for i in range(1, r + 1):  # filling the rest of the matrix
            for j in range(1, c + 1):
                f[i][j] = max(f[i - 1][j] + self.g, f[i][j - 1] + self.g,
                              f[i - 1][j - 1] + self.compare(a[i - 1], b[j - 1]))
        return f
```
-	EnumerateAlignments / EnumerateAlignments_lines (x, y, f, w, z): both are recursive methods that return a list of all possible alignments “ww” and “zz” of the “x” and “y” sequences. The slight difference between the two is that the first aligns two strings and the second aligns two lists of strings.
     - EnumerateAlignments (x, y, f, w, z): starts from the most right cell in the “f” matrix to the top most left cell, it back-traces from which previous cell was the current score cell derived (instead of making a tracing graph while creating “f” that will require space), and whenever a trace-back condition is fulfilled (the origin of the score in “f” is determined) the characters that correspond are concatenated to the sequences “w” and “z” (either a char or a “-”), and then sliced from the original strings; which are again fed back to the method, and so on till reaching the base case where both original strings become empty, then the resulting sequences are appended to the list of alignments. Since the origin of a score cell in “f” might be more than one, the recursion works on getting all possible alignments, because by the end of each case it goes back to the last step of recursion and checks for further matches and if a condition is fulfilled the whole process starts again.
         ```python
         def EnumerateAlignments(self, x, y, f, w, z):  # Gets a list of all alignments when x, y are strings
             i = len(x)
             j = len(y)
             if i == j == 0:  # The base case that ends the recursion
                 self.ww.append(w)
                 self.zz.append(z)
             if i > 0 and j > 0 and f[i][j] == f[i - 1][j - 1] + self.compare(x[i - 1], y[j - 1]):  # Diagonal result
                 self.EnumerateAlignments(x[0:i - 1], y[0:j - 1], f, x[i - 1] + w, y[j - 1] + z)
             if i > 0 and f[i][j] == f[i - 1][j] + self.g:  # Top result
                 self.EnumerateAlignments(x[0:i - 1], y, f, x[i - 1] + w, "-" + z)
             if j > 0 and f[i][j] == f[i][j - 1] + self.g:  # Left result
                 self.EnumerateAlignments(x, y[0:j - 1], f, "-" + w, y[j - 1] + z)
             return self.ww, self.zz
         ```
     - EnumerateAlignments_lines (x, y, f, w, z): works exactly like the first method but since the “w” and “z” are lists, it appends the strings that fulfill the condition to those lists (instead of simply concatenating them as characters).
     
         ```python
         def EnumerateAlignments_lines(self, x, y, f, w, z):  # Gets a list of all alignments when x, y are string lists
             i = len(x)
             j = len(y)
             if i == j == 0:  # The base case that ends the recursion
                 self.ww.append(w)
                 self.zz.append(z)
             if i > 0 and j > 0 and f[i][j] == f[i - 1][j - 1] + self.compare(x[i - 1], y[j - 1]):  # Diagonal result
                 self.EnumerateAlignments_lines(x[0:i - 1], y[0:j - 1], f, [x[i - 1]] + w, [y[j - 1]] + z)
             if i > 0 and f[i][j] == f[i - 1][j] + self.g:  # Top result
                 self.EnumerateAlignments_lines(x[0:i - 1], y, f, [x[i - 1]] + w, ["-"] + z)
             if j > 0 and f[i][j] == f[i][j - 1] + self.g:   # Left result
                 self.EnumerateAlignments_lines(x, y[0:j - 1], f, ["-"] + w, [y[j - 1]] + z)
             return self.ww, self.zz
         ```
Actually, the purpose of implementing the Needleman Wunsch algorithm though we’re working on the Hirschberg algorithm, is that it is one of base cases of the Hirschberg recursive algorithm.
   
 



