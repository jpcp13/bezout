To reproduce the experiment presented in the arXiv paper, you should have python and sage installed on your computer. Download and save this repository wherever you want on your computer and unzip the bezout-master.zip file; it should create a directory named bezout-master and containing two directories, one named python and one named tex.
Then follow the steps :

1. In file python/bezout_2.sage, check that :
    * the degree of the polynomial system is set to [2, 2, 2, 2] ; line 8 should write `deg = [2, 2, 2, 2]`
    * line 29 is uncommented ; this line should write 
        `P = load('P_'+''.join(str(e) for e in deg)+'.sobj')`
    * lines 30, 31 are commented, these lines should write 
        #~ P = [bz.rand_poly(n-1, m, deg, t, x) for i in range(n)] + xx
        #~ save(P, 'P_'+''.join(str(e) for e in deg))

2. From within the directory named python open a linux terminal;
    * from within this terminal execute the command `sage`; this will open a sage terminal
    * from within the sage terminal execute the command `runfile("bezout_2.sage")`; this will launch all the calculations; wait as long as necessary (the Groebner calculation takes about 2600s in this experiment)
    * After the calculations have been completed, check that the results are written in some text files with explicit names, located in the directory named tex/txt ; 


3. If you want these results to appear in the paper "Solving zero-dimensional polynomial systems: a practical method using Bezout matrices" ( file tex/bezout.pdf ), you must have latex installed on your computer; then
    * open a terminal in the directory named tex and execute the command `pdflatex bezout.tex`
    * this will create the file bezout.pdf, updated with the results of your own calculations.


