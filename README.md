# Additional exploration, code and ideas

This github repository is a collection of code and ideas related to the N-Queens problem which did not end up 
making it in my final thesis.

There is a lot more code written which is not in the thesis, but most of it is not very interesting, this github 
is for the potentially interesting results and code which did not make it into the final thesis.

Also, the plotting code and similar is relatively simple and easily reproducible and is therefore not included here or in the 
thesis.

## Alternative methods of computing $Q(n)$

A bunch of alternative methods to the standard backtracking method for computing $Q(n)$ are listed 
in the Qn folder. This includes an exponential runtime method $\sim O(36^n)$ and various optimizations on 
a method for calculating $Q(n)$ that I came up with that calculates $Q(n)$ by finding pairs of integers 
whose bitwise AND is zero, a GPU-implementation of this idea is also given.

Interestingly, finding the number of such pairs can be done in $o(|A||B|)$ time e.g. using k-d trees, where 
$A,B$ are the list of integers from which we draw pairs a in A and b in B and check if $a \& b = 0$. If 
the time-complexity can be significantly reduced from $O(|A||B|)$, computing e.g. $Q(28)$ would become quite feasible 
using the methods in the Qn folder.


## Other stuff
The other's folder contains various other interesting things such as 

An attempt at finding patterns in the $r$-symmetric n-queen solutions when varying $n$. 

A genetic algorithm. 

Attempts at finding a closed form solution.

3D rooks.

And more...



