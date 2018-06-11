# -*- coding: utf-8 -*-
"""
Created on Fri Jun  1 13:52:28 2018

@author: Ivan ALisson Cavalcante Nunes de Lima
"""

#!/Users/Rad/anaconda/bin/python
# (c) 2013 Ryan Boehning


'''A Python implementation of the Smith-Waterman algorithm for local alignment
of nucleotide sequences.

Ivan Alisson 2018:
    Modified class for alignment using Smith-Waterman and indexation of bases positions
'''
import Node as nd
import WorkFunctions as wFunc

class smithWaterman:
    # These scores are taken from Wikipedia.
    # en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm
    match    = 2
    mismatch = -1
    gap      = -1
    seq1     = None
    seq2     = None
    seq1Al   = None
    seq2Al   = None
    cut = None
    
    """
    mt = match score
    ms = miss match score
    gp = gap score

    sq1 = 1st sequence
    sq2 = 2nd sequence

    condition = True, show align.
                    False, hide align
                    
    writeAl = True, Write a file whith the alignment, for depuration purporses
            #### Set the path manually
    """
    def constructor(self, mt = 2, ms = -1, gp = -1, sq1 = None, sq2 = None, condition = False, writeAl=False):
        self.match    = mt
        self.mismatch = ms
        self.gap      = gp
        self.seq1     = sq1
        self.seq2     = sq2
        return self.main(condition, writeAl)


    def main(self, condition, writeAl):

        # The scoring matrix contains an extra row and column for the gap (-), hence
        # the +1 here.
        rows = len(self.seq1) + 1
        cols = len(self.seq2) + 1

        # Initialize the scoring matrix.
        score_matrix, start_pos = self.create_score_matrix(rows, cols)
        #print(start_pos, score_matrix)
        # Traceback. Find the optimal path through the scoring matrix. This path
        # corresponds to the optimal local sequence alignment.
        seq1_aligned, seq2_aligned, finalPosX, finalPosY = self.traceback(score_matrix, start_pos)

        assert len(seq1_aligned) == len(seq2_aligned), 'aligned strings are not the same size'
        # Pretty print the results. The printing follows the format of BLAST results
        # as closely as possible.
        alignment_str, idents, gaps, mismatches = self.alignment_string(seq1_aligned, seq2_aligned)
        alength = len(seq1_aligned)
        if condition:
            print()
            print(' Identities = {0}/{1} ({2:.1%}), Gaps = {3}/{4} ({5:.1%})'.format(idents,
                  alength, idents / alength, gaps, alength, gaps / alength))
            print()
            for i in range(0, alength, 60):
                seq1_slice = seq1_aligned[i:i+60]
                print('Seq01  {0:<4}  {1}  {2:<4}'.format(i + finalPosX, seq1_slice, i + len(seq1_slice)+finalPosX))
                print('             {0}'.format(alignment_str[i:i+60]))
                seq2_slice = seq2_aligned[i:i+60]
                print('Seq02  {0:<4}  {1}  {2:<4}'.format(i + finalPosY, seq2_slice, i + len(seq2_slice)+finalPosY))
                print()
        
        if writeAl:
            f = open("D:aligned.txt", "a")
            try:
                f.write("\n\nSeq01 {}\nStr00 {}\nSeq02 {} ".format(seq1_aligned, alignment_str, seq2_aligned))
            finally:
                f.close()

        self.cut = (idents / alength)*100
        self.seq2Al, self.seq1Al = wFunc.matchDegrees(self.seq2Al, self.seq1Al)
        
        return seq2_aligned, seq1_aligned, self.seq2Al, self.seq1Al

    def traceback(self, score_matrix, start_pos):
        '''Find the optimal path through the matrix.
        This 
        tion traces a path from the bottom-right to the top-left corner of
        the scoring matrix. Each move corresponds to a match, mismatch, or gap in one
        or both of the sequences being aligned. Moves are determined by the score of
        three adjacent squares: the upper square, the left square, and the diagonal
        upper-left square.
        WHAT EACH MOVE REPRESENTS
            diagonal: match/mismatch
            up:       gap in sequence 1
            left:     gap in sequence 2
        '''

        END, DIAG, UP, LEFT = range(4)
        aligned_seq1 = []
        aligned_seq2 = []
        seq1Al = []
        seq2Al = []
        x, y         = start_pos
        move         = self.next_move(score_matrix, x, y)
        while move != END:
            if move == DIAG:
                aligned_seq1.append(self.seq1[x - 1].getAmino())
                aligned_seq2.append(self.seq2[y - 1].getAmino())
                seq1Al.append(nd.node(self.seq1[x - 1].getAmino(), self.seq1[x - 1].getDegree()))
                seq2Al.append(nd.node(self.seq2[y - 1].getAmino(), self.seq2[y - 1].getDegree()))
                x -= 1
                y -= 1
            elif move == UP:
                aligned_seq1.append(self.seq1[x - 1].getAmino())
                aligned_seq2.append("-")
                seq1Al.append(nd.node(self.seq1[x - 1].getAmino(), self.seq1[x - 1].getDegree()))
                seq2Al.append(nd.node("-", -1))
                x -= 1
            else:
                aligned_seq1.append("-")
                aligned_seq2.append(self.seq2[y - 1].getAmino())
                seq1Al.append(nd.node("-", -1))
                seq2Al.append(nd.node(self.seq2[y - 1].getAmino(), self.seq2[y - 1].getDegree()))
                y -= 1
                

            move = self.next_move(score_matrix, x, y)
        
        aligned_seq1.append(self.seq1[x - 1].getAmino())
        aligned_seq2.append(self.seq2[y - 1].getAmino())
        seq1Al.append(nd.node(self.seq1[x - 1].getAmino(), self.seq1[x - 1].getDegree()))
        seq2Al.append(nd.node(self.seq2[y - 1].getAmino(), self.seq2[y - 1].getDegree()))
        
        self.seq1Al = list(reversed(seq1Al))
        self.seq2Al = list(reversed(seq2Al))
        return ''.join(reversed(aligned_seq1)), ''.join(reversed(aligned_seq2)), x-1, y-1

    def create_score_matrix(self, rows, cols):
        '''Create a matrix of scores representing trial alignments of the two sequences.
        Sequence alignment can be treated as a graph search problem. This function
        creates a graph (2D matrix) of scores, which are based on trial alignments
        of different base pairs. The path with the highest cummulative score is the
        best alignment.
        '''
        score_matrix = [[0 for col in range(cols)] for row in range(rows)]

        # Fill the scoring matrix.
        max_score = 0
        max_pos   = None    # The row and columbn of the highest score in matrix.
        for i in range(1, rows):
            for j in range(1, cols):
                score = self.calc_score(score_matrix, i, j)
                if score > max_score:
                    max_score = score
                    max_pos   = (i, j)

                score_matrix[i][j] = score

        assert max_pos is not None, 'the x, y position with the highest score was not found'

        return score_matrix, max_pos


    def calc_score(self, matrix, x, y):
        '''Calculate score for a given x, y position in the scoring matrix.
        The score is based on the up, left, and upper-left neighbors.
        '''
        similarity = self.match if self.seq1[x - 1].getAmino() == self.seq2[y - 1].getAmino() else self.mismatch

        diag_score = matrix[x - 1][y - 1] + similarity
        up_score   = matrix[x - 1][y] + self.gap
        left_score = matrix[x][y - 1] + self.gap

        return max(0, diag_score, up_score, left_score)


    def next_move(self, score_matrix, x, y):
        diag = score_matrix[x - 1][y - 1]
        up   = score_matrix[x - 1][y]
        left = score_matrix[x][y - 1]
        if diag >= up and diag >= left:     # Tie goes to the DIAG move.
            return 1 if diag != 0 else 0    # 1 signals a DIAG move. 0 signals the end.
        elif up > diag and up >= left:      # Tie goes to UP move.
            return 2 if up != 0 else 0      # UP move or end.
        elif left > diag and left > up:
            return 3 if left != 0 else 0    # LEFT move or end.
        else:
            # Execution should not reach here.
            raise ValueError('invalid move during traceback')


    def alignment_string(self, aligned_seq1, aligned_seq2):
        '''Construct a special string showing identities, gaps, and mismatches.
        This string is printed between the two aligned sequences and shows the
        identities (|), gaps (-), and mismatches (:). As the string is constructed,
        it also counts number of identities, gaps, and mismatches and returns the
        counts along with the alignment string.
        AAGGATGCCTCAAATCGATCT-TTTTCTTGG-
        ::||::::::||:|::::::: |:  :||:|   <-- alignment string
        CTGGTACTTGCAGAGAAGGGGGTA--ATTTGG
        '''
        # Build the string as a list of characters to avoid costly string
        # concatenation.
        idents, gaps, mismatches = 0, 0, 0
        alignment_string = []
        for base1, base2 in zip(aligned_seq1, aligned_seq2):
            if base1 == base2:
                alignment_string.append('|')
                idents += 1
            elif '-' in (base1, base2):
                alignment_string.append(' ')
                gaps += 1
            else:
                alignment_string.append(':')
                mismatches += 1

        return ''.join(alignment_string), idents, gaps, mismatches

    def getCut(self):
        return self.cut