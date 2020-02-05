# CS362-HW3
  
1. Our code is python, so no need to compile, and has no dependencies besides
   sys

2. There are currently no known bugs in the code.

3. The order of tie-breaking for alignments matters a bit; whether, all else
   equal, the algorithm chooses to move horizontally, diagonally, or
   vertically in the final alignments. Our order is horizontal movements
   are prioritized, followed by diagonal, then vertical; we chose this to match
   the example output.
   We also included a helper function in both programs entitled PrintNice()
   which takes the algorithm output and prints it according to the
   specifications set by the example output.

4. This project was written by Daniel Busis and Kristin Albright. We consulted
   with Nikko Baer to confirm we were getting the same results on the same
   input.
