import nmw


class hb:

    def __init__(self, g, m, d, t):
        self.g = g
        self.d = d
        self.m = m
        self.t = t

    def compare(self, x, y):  # Compares sequence elements x and y to determine the match/difference score
        if x == y:
            return self.m
        else:
            return self.d

    def ComputeAlignmentScore(self, x, y):
        L = [j * self.g for j in range(len(y) + 1)]  # Filling the first row of the Alignment Score matrix
        K = [0 for _ in range(len(y) + 1)]  # K as the second row
        for i in range(1, len(x) + 1):
            L, K = K, L  # Now L will be the current row and K is the previous line (Previously calculated)
            L[0] = i * self.g  # the first cell is always calculated j * g
            for j in range(1, len(y) + 1):  # Calculating the current line from top and diagonal previous cells and left cell from L
                L[j] = max(L[j - 1] + self.g, K[j] + self.g,
                           K[j - 1] + self.compare(x[i - 1], y[j - 1]))  # Based on the principle of optimality
        return L

    def Hirschberg(self, x, y):  # when the sequences are characters
        if len(x) == 0:  # If x is empty, y is aligned with "-"
            WW = ['-' * len(y)]
            ZZ = [y]
        elif len(y) == 0:  # If y is empty, x is aligned with "-"
            WW = [x]
            ZZ = ['-' * len(x)]
        elif len(x) == 1 or len(
                y) == 1:  # If one of their lengths is 1, Needleman does the job with a matrix: 1 x len(y), or: len(x) x 1
            align_nmw = nmw.nmw(self.g, self.m, self.d)
            f = align_nmw.F(x, y)
            z = ""
            w = ""
            WW, ZZ = align_nmw.EnumerateAlignments(x, y, f, w, z)
        else:
            i = len(x) // 2
            Sl = self.ComputeAlignmentScore(x[0:i], y)  # Computing the AS of the first half of x with y
            Sr = self.ComputeAlignmentScore(x[i:len(x)][::-1], y[::-1])  # Computing the AS of the reversed second half of x with  reversed y
            # Due to the Reversion we get exactly the last line of the first half and the first line of the second half
            S = [k + t for k, t in zip(Sl, Sr[::-1])]  # Sum of the Sl and reversed Sr
            m = max(S)  # Maximum value of S to find where to slice y for the next alignment
            J = [i for i, j in enumerate(S) if j == m]  # Getting the indices of the m occurrence in S
            WW = []
            ZZ = []
            for j in J:  # To slice y with all possible solutions
                if self.t:
                    print(i, ",", j)
                WWl, ZZl = self.Hirschberg(x[0:i], y[0:j])
                WWr, ZZr = self.Hirschberg(x[i:len(x)], y[j:len(y)])
                for z, c in zip(WWl, ZZl):  # Just to make all possible permutations of the outcome
                    for k, e in zip(WWr, ZZr):
                        WW.append(z + k)
                        ZZ.append(c + e)
        return WW, ZZ

    def Hirschberg_lines(self, x, y):  # when sequences are lines, same as first one but some changes to deal with it as a list
        if len(x) == 0:  # If x is empty, y is aligned with "-"
            print("x", x)
            WW = ['-' * len([y])]
            ZZ = [y]
        elif len(y) == 0:  # If y is empty, x is aligned with "-"
            print("y", y)
            WW = [x]
            ZZ = ['-' * len([x])]
        elif len([x]) == 1 or len([y]) == 1:  # If one of their lengths is 1, Needleman does the job with a matrix: 1 x len(y), or: len(x) x 1
            align_nmw = nmw.nmw(self.g, self.m, self.d)
            f = align_nmw.F(x, y)
            z = []
            w = []
            WW, ZZ = align_nmw.EnumerateAlignments_lines(x, y, f, w, z)
        else:
            i = len(x) // 2
            Sl = self.ComputeAlignmentScore(x[0:i], y)  # Computing the AS of the first half of x with y
            Sr = self.ComputeAlignmentScore(x[i:len(x)][::-1], y[::-1])  # Computing the AS of the reversed second half of x with  reversed y
            # Due to the Reversion we get exactly the last line of the first half and the first line of the second half
            S = [k + t for k, t in zip(Sl, Sr[::-1])]  # Sum of the Sl and reversed Sr
            m = max(S)  # Maximum value of S to find where to slice y for the next alignment
            J = [i for i, j in enumerate(S) if j == m]  # Getting all the indices of the m occurrence in S
            WW = []
            ZZ = []
            for j in J:  # To slice y with all possible solutions
                if self.t:
                    print(i, ",", j)
                WWl, ZZl = self.Hirschberg_lines(x[0:i], y[0:j])  # left part
                WWr, ZZr = self.Hirschberg_lines(x[i:len(x)], y[j:len(y)])  # right part
                print(WWl, WWr, "\n", ZZl, ZZr)
                for z, c in zip(WWl, ZZl):  # Just to make all possible permutations of the outcome
                    for k, e in zip(WWr, ZZr):
                        WW.append([z].extend([k]))  # extend instead of simple concatenation
                        ZZ.append([c].extend([e]))
        return WW, ZZ
