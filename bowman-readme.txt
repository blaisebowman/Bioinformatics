Blaise Bowman, CIS 4930: Bioinformatics Project #1

To run this program, I transferred bowman.cpp and sequenceGenerator.cpp to Linux via WinSCP. 

Within Linux, I used the following commands, in order:

g++ sequenceGenerator.cpp 
./a.out
	//The above generates the pairs of sequences with a given length, which are displayed in files 50characters.o1 -> 5000characters.o8

g++ bowman.cpp 
	//Compile the program for execution

./a.out assignment-1.fasta //replace "assignment-1.fasta" with the file for testing 

//Algorithm Analysis Below:
Note: to test algorithm analysis, you must comment out Part F's code as explained in bowman.cpp

time ./a.out 50characters.o1
	//The above command determines the execution time of the algorithm with the sequence pair 50 characters in length is passed 
valgrind --tool=massif ./a.out 50characters.o1 
	//The above command determines the total memory usage when executing the program, find memory usage at the peak
ms_print massif.out.1558

time ./a.out 100characters.o2
valgrind --tool=massif --pages-as-heap=yes ./a.out 100characters.o2
ms_print massif.out.1571

time ./a.out 250characters.o3
valgrind --tool=massif ./a.out 50characters.o3 
ms_print massif.out.1577


time ./a.out 500characters.o4
valgrind --tool=massif ./a.out 500characters.o4
ms_print massif.out.1583


time ./a.out 1000characters.o5
valgrind --tool=massif ./a.out 1000characters.o5
ms_print massif.out.1592


time ./a.out 1500characters.o6
valgrind --tool=massif ./a.out 1500characters.o6 
ms_print massif.out.1595


time ./a.out 2000characters.o7
valgrind --tool=massif ./a.out 2000characters.o7
ms_print massif.out.1598


time ./a.out 10000characters.o8
valgrind ./a.out 10000characters.o8
ms_print massif.out.1669


Note: Testing for part F works when determining the best alignment and if there are multiple alignments. For longer DNA sequences (50 characters or more),
my computer freezes up when actually determining all possible optimal alignments. Though, I believe this was out of the scope of the project. 
