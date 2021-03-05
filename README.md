# Dynamic Programming and Pairwise Sequence Alignment

The program performs pairwise sequence alignment of two protein sequences or two DNA sequences depending on the user's configuration of parameters. It can also perform global sequence alignment,  semi-global sequence alignment, and local sequence alignment depending on the user's configuration. The program chooses the corresponding scoring matrix - dnaMatrix or BLOSUM45 - file depending on if it is a protein sequence alignment or DNA sequence alignment.

  Parameters for the program:
    -i First sequence [File In and should be FASTA format]
    -j Second sequence [File In and should be FASTA format]
    -p Align protein sequences or DNA sequences [T/F]
    -atype global, or semi-global [G/S]
    -o alignment output file [File Out]

  Possible ways of running the program - align.py. 
  
    • $python3 align.py -i seq1.txt -j seq2.txt -o out.txt -p F -atype G
      Perform global pairwise sequence alignment of two DNA sequences in seq1.txt and
      seq2.txt. The output will be in out.txt.

    • $python3 align.py -i seq1.txt -j seq2.txt -o out.txt -p T -atype S
      Perform semi-global pairwise sequence alignment of two protein sequences in seq1.txt and
      seq2.txt. The output will be in out.txt.

  The order of the parameters does not matter. For example:
    
    $python3 align.py -i seq1.txt -j seq2.txt -o out.txt -p F -atype G
    &
    $python3 align.py -i seq1.txt -j seq2.txt -p F -atype G -o out.txt
