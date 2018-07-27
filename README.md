If you want to reproduce the experiment presented in the arXiv paper, you should have python and sage installed on your computer. Download and save this repository at a place of your choice on your computer and unzip the bezout-master.zip file; it should create a directory named bezout-master and containing two directories, one named python and one named tex.
Then follow the steps :

1. In the file python/bezout_2.sage :
    * do not modify line 8 ; this line should write `deg = [2, 2, 2, 2]`
    * uncomment line 30, this line should write `P = load('P_'+''.join(str(e) for e in deg)+'.sobj')`
    * comment line 31, this line should write `# P = [bz.rand_poly(n-1, m, deg, t, x) for i in range(n)] + xx`

2. From within the directory named python open a linux terminal;
    * from within this terminal execute the command `sage`; this will open a sage terminal
    * from within the sage terminal execute the command `runfile("bezout_2.sage")`; this will make all the calculations; wait as long as necessary.
    * After the calculations have been completed, check that the results are written in some text files with explicit names, located in the directory named tex/txt ; 


3. If you want these results to appear in the paper "Solving zero-dimensional polynomial systems: a practical method using Bezout matrices" ( file tex/bezout.pdf ), you must have latex installed on your computer; then
    * open a terminal in the directory named tex and execute the command `pdflatex bezout.tex`
    * this will create the file bezout.pdf, updated with the results of your own calculations.


