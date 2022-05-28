# Code source Structure 
The Program consists of 3 parts: the Main **(Hirschberg.py)**, and 2 modules: the hb  **(Hirschberg) class (hb.py)**, and the **nmw (Needleman-Wunsch) class (nmw.py)**

- **The hb class:** For the creation of “hb” objects, after which we call its Hirschberg method  to create the list of aligned sequences. 
- **The nmw class:** For the creation of “nmw” objects, after which we call its  EnumerateAlignments method to create the list of aligned sequences.
- **The main:** The starting point of the program’s execution. Results in both the Needleman Wunsch and Hirschberg Alignment.

## The Main (Hirschberg.py): 
**The program takes 8 arguments:**

 1.03 optional arguments: 
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
   
 



