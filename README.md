## Motivation

The basic idea in this implementation is to use the Lenstra-Lenstra-Lovasz ([LLL](https://en.wikipedia.org/wiki/Lenstra%E2%80%93Lenstra%E2%80%93Lov%C3%A1sz_lattice_basis_reduction_algorithm)) algorithm to compute an approximate linear dependence (with integer coefficients) between powers of a real root (obtained via Newton's numerical method) by means of a lattice reduction procedure.


### Code summary

0. Code is divided up into 6 sections: 
	1. Basic arithmetic operations for lattices; 
	2. Class of Rank 2 Lattices; 
	3. Class of high rank lattices; 
	4. Basic arithmetic of polynomials with integer coefficients; 
	5. a version of Newton's method; 
	6. the Main method with basic UI

1. If you save then import the file into Python, you may run the main() method by passing as argument the polynomial you wish to factor written as an array of coefficients:
	[a_0, a_1, ..., a_k]
where a_0 will be the constant term and a_k the coefficient of the highest degree term.

2. If you execute without argument to *main()*, you will be prompted to enter the degree and coefficients of polynomials until you abort the process. 

3. This algorithm relies on an initial pass with Newton's method, which has not been numerically optimized to the highest possible standard. As such, the factorization works best if all the roots are real.

4. If input polynomial has complex roots, the algorithm runs correctly but with the currently enforced limit on the number of iterations allowed, the algorithm may not yield a full factorization when the coefficients and/or the degree are not very small. 

5. Modifications can be made to optimize the code by using a stable Newton method which works equally well over complex numbers (with arbitrary precision) and by increasing the number of allowed iterations depending on the size of the input. 

6. Some examples you can try: [1, -2, 1], [1, 4, 6, 4, 1], [-1, 0,0,0, 1], [1, 0, 2, 0, 1] or larger examples like:
	* [4147, 6090, -116, -172, -210, 4, 1], factors into (x^3 - 29)(x^2+17x+11)(x-13),
	* [-5436640836, 95024916, -447033, 152, 1] as (x-123)(x-438)(x+834)(x-121)
	* [156, -9, -230, -2, 1] as (x^2 + 15x + 12)(x^2 - 17x + 13)

Alternatively, you may generate test polynomials to factor with MultiplyPolys(p1, p2) implemented in section 4 of the file.
