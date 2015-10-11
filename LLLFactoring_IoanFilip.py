### IOAN FILIP, COLUMBIA UNIVERSITY, MATH PHD

### 7/11/15

### RECURSE CENTER WINTER 2, 2016 BATCH APPLICATION 

######################################################################

### Polynomial Factoring up to Q-units using LLL Lattice Reduction ###

######################################################################

### 0. Some useful Imports


import math, cmath
from random import random
from fractions import Fraction, gcd

### Some global Variables to fix error bound and iterations allowed for Newthon's method

gError = 0.000000000001
MaxNewtonIterates = 100
NewtonTrials = 50


######################################################################

### 1. Basic Arithmetic Operations for Lattices


def nearestInt(x):
	m = 0
	
	if (abs(x - int(x)) >= 0.5):
		if x>=0:
			m = int(x) + 1
		else:
			m = int(x) - 1
	else:
		m = int(x)

	return m


def normQ(N, alpha_vector, b):
	inner_sum = 0
	outer_sum = 0
	for indx in range(0, len(alpha_vector)):
		inner_sum += alpha_vector[indx]*b[indx]
		outer_sum += abs(b[indx])**2

	norm_result = outer_sum + N*abs(inner_sum)**2
	return norm_result


def add(b1, b2):
	sum = []
	for i in range(0, len(b1)):
		sum.append(b1[i] + b2[i])

	return sum

def innerQ(N, alpha, b1, b2):
        sum_vector = add(b1, b2)
        difference_norms = normQ(N, alpha, sum_vector) - normQ(N, alpha, b1) - normQ(N, alpha, b2)
        return 0.5*difference_norms


def scale(a, b):
	scaled = []
	for i in range(0, len(b)):
		scaled.append(a*b[i])
	return scaled 

def GramSchmidt(vectors, N, alpha):
	rank = len(vectors)
	vectors_new = []
	
	vectors_new.append(vectors[0])

	for i in range(1, rank):
		new_vector = vectors[i]

		for j in range(0, i):
			c_ij = innerQ(N, alpha, new_vector, vectors_new[j])/normQ(N, alpha, vectors_new[j])
			new_vector = add(new_vector, scale(-c_ij, vectors_new[j]))

		vectors_new.append(new_vector)

	return vectors_new


def Succ_Dist(vectors, N, alpha):
	# Return successive lengths of 'vectors' realized as flag, namely disc(Li/Li-1), namely length of proj of vi on sublattice Li-1

	dist = []
	new_vectors = GramSchmidt(vectors, N, alpha)
	for i in range(0, len(vectors)):
		dist.append(math.sqrt(normQ(N, alpha, new_vectors[i])))

	return dist

########################################################################

### 2. Implementation of Rank 2 Lattice Class -- Base case and recursive step



class Rk2_Lattice(object):
		
	def __init__(self, v1, v2, alpha, N):
		self.N = N
		self.alpha = alpha
		
		norm_v1 = normQ(N, alpha, v1)
		norm_v2 = normQ(N, alpha, v2)
		
		if norm_v2 >= norm_v1:
			self.b1 = v1
			self.b2 = v2
		else:
			self.b1 = v2
			self.b2 = v1
			
		inn_b1b2 = 2*innerQ(N, alpha, self.b1, self.b2) 
		# This is 2<b1, b2> = -2*N*alpha
		
		self.reduced = abs(inn_b1b2) <= norm_v1 and norm_v1 <= norm_v2
		
		
	def has_basis(self):
		inn_b1b2 = 2*innerQ(self.N, self.alpha, self.b1, self.b2)
		return not (inn_b1b2**2 - 4*normQ(self.N, self.alpha, self.b1)*normQ(self.N, self.alpha, self.b2) == 0)
	
	def is_reduced(self):
		return self.reduced
 	
			
	def discriminant(self):
		inn_b1b2 = 2*innerQ(self.N, self.alpha, self.b1, self.b2)
		norm_v1 = normQ(self.N, self.alpha, self.b1)
		norm_v2 = normQ(self.N, self.alpha, self.b2)
		return math.sqrt(abs(norm_v1*norm_v2 - (inn_b1b2**2)/4.0))


	def reduce(self):
		if self.is_reduced():
			return None		
		elif normQ(self.N, self.alpha, self.b1) > normQ(self.N, self.alpha, self.b2):
			v1 = self.b1
			v2 = self.b2

			self.b1 = v2
			self.b2 = v1

		else:
			inn_b1b2 = 2*innerQ(self.N, self.alpha, self.b1, self.b2)
			n = inn_b1b2/(2*normQ(self.N, self.alpha, self.b1))
			
			m = nearestInt(n)
						
			self.b2 = add(self.b2, scale(-m, self.b1))
			
			
			# After this calculation, SHOULD already have that abs(self.inn_b1b2) is at most q(b1); condition below checks the reduced property of new basis
			if(normQ(self.N, self.alpha, self.b2) >= normQ(self.N, self.alpha, self.b1)):
				self.reduced = True
			# Sets reduced variable to True; any discrepancy with the inequalities reduced check results from calculations exceeding error bounds 
				return None

		self.reduce()		

	def return_basis(self):
		return [self.b1, self.b2]


##############################################################################


### 3. Implementation of Rank n Lattice Class for degree n polynomial factorization
	    

class Rkn_Flag(object):

        ### Sets an absolute upper bound on the number of iterations of the basis reduction algorithm
        ### For debugging purposes only: expect number is polynomial in the log-size of the flag
        
        reductions_max = 50

        def __init__(self, vectors, N, alpha):
                # Creates Flag from list of vectors, ASSUMING the associated sublattices are pure (if not so, any torsion part is ignored in all calculations -- possible source of ERROR?)
                # Also assumes length of vectors equals that of alpha equals rank of flag

                self.N = N
                self.alpha = alpha
                self.vectors = vectors
                self.rank = len(vectors)

        def get_distances(self):
                return Succ_Dist(self.vectors, self.N, self.alpha)

        def disc(self):

                # Computes the discriminant of the Lattice determined by the full set of basis vectors

                succ_dist = self.get_distances()

                disc = 1

                for i in range(0, self.rank):
                        disc = disc*succ_dist[i]

                return disc

        def size(self):

                succ_dist = self.get_distances()

                size = 1
                for i in range(0, self.rank):
                        size = size*succ_dist[i]**(self.rank - i)

                return size

        def is_cReduced(self, c):
                # Checks if Flag is c-reduced for a parameter c and returns the first index of the subflag which isn't
                index_list = []

                c_reduced = False

                succ_dist = self.get_distances()

                for i in range(0, self.rank - 1):
                        if succ_dist[i+1]**2 < (succ_dist[i]**2)/c:
                                index_list.append(i)

                if(len(index_list) == 0):
                        return [True, index_list]
                else:
                        return [False, index_list]

        
        def size_Reduce(self):
                # Similar to Gram-Schmidt procedure where just replace with nearest int, not projection coeff itself
                # Every vector in Size-Reduced flag satisfies the property that the projection to the smaller rank sublattice lies in the latter's fundamental domain

                rank = len(self.vectors)

                vectors_new = []

                vectors_new.append(self.vectors[0])

                for i in range(1, rank):
                        new_vector = self.vectors[i]

                        for j in range(0, i):
                                c_ij = nearestInt(innerQ(self.N, self.alpha, self.vectors[i], vectors_new[j])/normQ(self.N, self.alpha, vectors_new[j]))
                                new_vector = add(new_vector, scale(-c_ij, vectors_new[j]))

                        vectors_new.append(new_vector)

                self.vectors = vectors_new

        def Basis_Reduce(self, c):		
		# Flag reduction just iterates rank 2 reduction of consecutive quotients by pure sublattice (no torsion) with same normQ-with parameters
		# Assumes reduction algorithm is called with parameter c = 2, to ensure polynomial time and compatibility with rank 2 iteration 	
                
                reductions_count = 0
		
		# There are only finitely many flags of smaller size than given self, so the algorithm will halt: every iteration strictly reduces the size of flag
                while(not self.is_cReduced(c)[0] and reductions_count < Rkn_Flag.reductions_max):
                        index_list = self.is_cReduced(c)[1]
			
                        self.size_Reduce()				
                        
                        while(len(index_list) != 0):
                                index = index_list[0]
                                Rank2_quotient = Rk2_Lattice(self.vectors[index], self.vectors[index+1], self.alpha, self.N)
                                if Rank2_quotient.is_reduced():
                                        index_list.pop()
                                else:
                                        Rank2_quotient.reduce()
                                        break

                        if len(index_list) == 0:
                                break
			
                        new_quotient = Rank2_quotient.return_basis()

                        self.vectors[index] = new_quotient[0]
                        self.vectors[index+1] = new_quotient[1]

                        reductions_count += 1
                        

        def short_vector(self):
		# Applyies the Flag basis-reduction algorithm with c = 2, which is the standard choice
                self.Basis_Reduce(2)
		
		# A shortest vector will be the first in a c=2-Reduced Flag	
                return self.vectors[0]


#################################################################################


### 4. Implementation of polynomials with rational coefficients as arrays
###    plus Implementation of basic polynomial arithmetic operations


def deg(polynomial):
    global gError
    
    d = len(polynomial)

    while d > 0:
        d -= 1
        if(math.fabs(polynomial[d]) > gError):
            break

    return d
    

def truncate(poly, degree):    
    return poly[0:degree + 1]

def printPoly(poly):
    global gError
    poly_to_str = ""

    for i in range(len(poly)):
        if(math.fabs(poly[i] - 0) < gError):
            continue
        else:
            poly_to_str += str(poly[i])+" * x^"+str(i)+" + "

    return poly_to_str[0: len(poly_to_str) - 3]

def isZero(poly):
    global gError
    
    for i in range(len(poly)):
        if math.fabs(poly[i]) > gError:
            return False

    return True


# Division of Integer coefficient poly over Q: returns q and r such that f = q*g + r, deg(r)<deg(g) up to Q-units
# Implemented Recursively

def DivWithRemainder(f, g):
    global gError
    global deg

    if isZero(g):
        raise ValueError("Division by zero polynomial")

    if isZero(f):
        return [[], [0]]
    
    d1 = deg(f)
    d2 = deg(g)

    
    if d1 < d2:
        return [[0], f]

    if d2 == 0:
        poly_diff = []
        for i in range(d1+1):
            poly_diff.append(f[i]/g[0])
        return [poly_diff, [0]]

    diff = d1 - d2

    poly_diff = []
    coeff = f[d1]/g[d2]

    for i in range(diff):
        poly_diff.append(0)
        
    for i in range(d2 + 1):
        poly_diff.append(g[i]*coeff)

    subtract_polys = []

    for i in range(d1 + 1):
        subtract_polys.append(f[i] - poly_diff[i])

    subtract_polys = truncate(subtract_polys, deg(subtract_polys))

    
    [inductive_poly, remainder] = DivWithRemainder(subtract_polys, g)

    
    deg_difference = deg(inductive_poly)

    zero_shift = 0

    if len(inductive_poly)== 0:
        zero_shift = diff
    else:
        zero_shift = diff - deg_difference - 1

    for i in range(zero_shift):
        inductive_poly.append(0)

    inductive_poly.append(coeff)
    
    
    return [inductive_poly, remainder]
    
    
def MultiplyPolys(poly1, poly2):
    global deg
    global gError
    
    d1 = deg(poly1)
    d2 = deg(poly2)

    mult = []
    
    for i in range(d1+d2+1):
        mult.append(0)

    for i in range(d1 + 1):
        for j in range(d2 + 1):
            mult[i+j] += poly1[i]*poly2[j]

    mult_degree = deg(mult)

    return truncate(mult, mult_degree)

    
### Function to deal with Polynomials and floating point coefficients turning them into integer coefficients UP to units of Q,
# i.e. polynomial expression normalized over Q so that coefficients are integral


def IntegralCoeff(poly):
    global gError
    global deg
    
    non_zeros = []

    for i in range(len(poly)):
        if poly[i] > gError:
            non_zeros.append(poly[i])

    denominators = []
    product = 1
    
    for i in range(len(non_zeros)):
        poly[i] = Fraction(poly[i]).limit_denominator()
        if poly[i].denominator != 1:
            denominators.append(poly[i].denominator)
            product *= poly[i].denominator

    lcm = 1
    
    if len(denominators) == 0:
        lcm  = product
    else:
        lcm = int(product/GCD(denominators))

    for i in range(len(poly)):
        temp_fraction = Fraction(poly[i]*lcm).limit_denominator()
        poly[i] = temp_fraction.numerator

    # Return vector keeps track of the normalization fraction to scale the polynomial
    return [poly, Fraction('1/'+str(lcm))]


def GCD(array):
    k = len(array)
    if k == 0:
        return 0
    elif k == 1:
        return array[0]
    else:
        return gcd(GCD(array[0:int(k/2)]), GCD(array[int(k/2):k]))

def Evaluate(poly, float_num):
    global deg, gError
    
    value = 0

    for k in range(deg(poly)+1):
        value += poly[k]*float_num**k

    return value


def GetFactor(poly):
    global deg, gError

    degree = deg(poly)
    
    absolute_coeffs = []

    for i in range(degree + 1):
            absolute_coeffs.append(math.fabs(poly[i]))

    ### THIS IS A CRUCIAL ISSUE: normalization of parameters N and epsilon to work at every scale & return correct integral factors

    M = max(absolute_coeffs) + 2

    N1 = ((16/9)*M)**(degree) + 105
    N2 = 11**10 + (3*M)**2 

    N = min(N1, N2)
    ## TO DO: try other cases depending on growth of coefficients and degrees

    epsilon = 1/(N**(1.1) + 1)

    [alpha, Newton_result] = Newton(poly, min(gError, epsilon))

    if not Newton_result:
        return [None, poly]

    alphas = [1]

    multiplier = 1

    for i in range(degree-1):
        alphas.append(alpha*multiplier)
        multiplier *= alpha

    # Initiates the standard basis for rank = degree flag
    vectors = []

    for i in range(degree):
        new_vec = []
        for j in range(degree):
            new_vec.append(0)
        new_vec[i] = 1

        vectors.append(new_vec)

    # Initiates Flag with given data

    Poly_lattice = Rkn_Flag(vectors, N, alphas)

    # Returns shortest vector in the lattice, interpreted as a POLYNOMIAL:

    shortest_solution = Poly_lattice.short_vector()

    #print("Iteration has reached a short vector solution of Rk n lattice: ")
    #print(shortest_solution)

    ## NEXT: check if the shortest vector indeed provides a factor of the polynomial:

    return [False, shortest_solution]

    
#########################################################################


### 5. Implementation of Newton's method for root finding over Complex Numbers
###     This is needed to initiate the Lattice reduction problem 


def Newton(poly, epsilon):
    global deg, gError, MaxNewtonIterates, NewtonTrials
# Takes in a polynomial of degree d and an error; returns root within error
# Polynomial given as a d+1 array of floats by F(x) = a_{0}x^0 + ... + a_{d}x^d
    
    d = deg(poly)
    
    constant = (d == 0)
        
    if(constant and (poly[0] == 0)):
        return [0, True]
    elif(constant):
        #print("Constant non-zero Polynomial!")
        return [1, False]
 
    x_n = complex(random(), random())
    x_m = complex(random(), random())

    # In case Newtown doesn't converge to solution from inital value of -1
    i = 0
    Tr = 0
    while (i<MaxNewtonIterates and Tr<NewtonTrials):
        poly_at_x_m = 0
        derivative_at_x_m = 0

        for j in range(0,d+1):
            poly_at_x_m = poly_at_x_m + poly[j]*(x_m)**j
            if(j == 0):
                continue
            else:
                derivative_at_x_m = derivative_at_x_m + j*poly[j]*(x_m)**(j-1)

        if(abs(derivative_at_x_m) < epsilon):
            if(abs(poly_at_x_m) < epsilon):
                return[x_m, True]
# Slightly perturb solution so as not to divide by zero
            x_m = x_m + 0.0001*complex(random(), random())
            i = i + 1
            continue

        x_n = x_m
        x_m = x_m - poly_at_x_m/derivative_at_x_m

        if abs(Evaluate(poly, x_m)) < epsilon:
            break

        if i < MaxNewtonIterates - 2:
            i = i+1
        else:
            x_n = complex(random(), random())
            x_m = complex(random(), random())
            i = 0
            Tr += 1
            

    if(abs(Evaluate(poly, x_m)) < epsilon):
        return [x_m, True]
    else:
        return [x_m, False]


###########################################################################


### 6. MAIN: basic UI for inputing a polynomial, attempts factorization via LLL and outputs result
###     It has following features:
###     A Flag construction given polynomial
###     Factorization algorithm,applied recursively, while degree > 1


###import sys


def main(poly = None):
    global deg
    global gError
    
    start_input = False
    
    if poly is None:
        user_in = True
    else:
        user_in = False
        start_input = True

    while(user_in or start_input):
        if not start_input:
            print("Enter integer coefficient polynomial to factorize: ")
            degree = int(input("What degree? "))
            poly = []
        
            for indexx in range(degree + 1):
                poly.append(int(input("Next coefficient... x^"+str(indexx)+" : ")))

        degree = deg(poly)
        
        poly = truncate(poly, degree)

        print("You have entered the polynomial: "+printPoly(poly))

        print("Performing factorization over Q... ")

        ######## Include full LLL-based algo HERE: appeal to recursive getFactor function ###########

        if degree > 0:
            [irreducible, next_factor] = GetFactor(poly)
            if irreducible is None:
                return "Newton's method failed to find a good initial guess. Try again..."
        else:
            next_factor = poly

        factors = []
        scales = []

        if deg(next_factor) == degree or deg(next_factor) == 0:
            irreducible = True

        polyA = poly
        
        while(not irreducible):
            integral_coeff_temp = IntegralCoeff(next_factor) 
            factors.append(integral_coeff_temp[0])
            scales.append(integral_coeff_temp[1])
            rest_poly_fact_temp = polyA
            [polyA, B] = DivWithRemainder(polyA, next_factor)
            polyA = truncate(polyA, deg(polyA))
            ### CHECK about the remainders at every step: consistence, should be zero!!
            if not isZero(B):
                irreducible = True
            # Would report a faulty factorization at this stage; the last 2 factors in the final answer are WRONG if this ever occurs

                factors.pop()
                scales.pop()
                polyA = rest_poly_fact_temp
                break

            if deg(polyA) == 1 or deg(polyA) == 0:
                irreducible = True
                continue

            
            [irreducible, next_factor] = GetFactor(polyA)

            if deg(next_factor) == deg(polyA) or deg(next_factor) == 0:
                irreducible = True

        if deg(polyA) < degree or len(factors) == 0:
            factors.append(IntegralCoeff(polyA)[0])
            scales.append(IntegralCoeff(polyA)[1])

        ###############################################################################################

        print("Your polynomial was "+printPoly(poly))
        if len(factors) == 1:
            print("This polynomial is irreducible over Q.")
        else:
            print("The factors are: ")

            for i in range(len(factors)):
                if deg(factors[i]) > 0:
                    print(printPoly(factors[i]))

            print("Checking that the product of the factors give initial poly, up to Q-units: ")

            product_poly = [1]
            product_scales = 1

            for i in range(len(factors)):
                # Factors already integer coefficient; no need for IntegralCoeff function
                product_poly = MultiplyPolys(product_poly, factors[i])
                product_scales = (product_scales*scales[i]).limit_denominator()

            print("The product of the factors is: "+printPoly(product_poly)+" with a rational scaling by "+str(product_scales))           


        ### User input to continue Factoring Polynomials
        if user_in:
            user = input("Continue?... y/n ")
            if user == "n" or user == "N":
                user_in = False
        else:
        ### Call at input control
            start_input = False



if __name__ == "__main__":
    main()
