for specific information please see http://www.bennettlovelady.net/software/finding-a-better-truss/

solver.py is where the juicy bits lie

there are 7 main steps within the solver program:

step 1: read the input data
step 2: calculate the sin and cos coefficients for each member
step 3: construct the matrix that represents the system
step 4: solve the system of equations to find the force in each member
step 5: find the maximum allowed load
step 6: output the results


--- step 1 ---
the input data are stored in:
	"joints.in": contains the x,y coordinates of each joint in mm
	"members.in": says which two joints each member connects to
	"applied.in": is the list of applied forces at which joints
	"supports.in": is the two end points, used for calculating reactions
	"material.in": contains the material properties - E, yield stress etc


--- step 2 ---
in any equation that a member "k" contributes to, the equation appears as:
... + cos(theta) . F(k) + ... = 0
or sin(theta).
This step calculates the cos() and sin() coefficients for each member.


--- step 3 ---
the truss is a system of linear equations, that is all the forces are multiplied by a constant value cos(theta) or sin(theta). 
This can be represented as a matrix equation C.Q=P, where:
	C is a 2nx2n matrix of coefficients, n being the number of joints
	Q is a vector of unknowns with 2n elements = (2n-3) members + 3 reactions
	P is a vector of applied or known forces

The matrix C can be made to represent the system of equations by placing the
cos() and sin() values from before in the right positions of the matrix


--- step 4 ---
The system is solved by inverting the matrix C:
	C . Q = P
	so Q = C^-1 . P

In the case where the matrix C is singular (has no inverse), this can be solved by numerical approximation by least squares analysis.

The result of this is that Q contains all our unknown forces


--- step 5 ---
In this step we go through each member and calculate the maximum load before that member breaks

1) yield.
To test yield we have to make sure that:
	F(k) / A <= stress_allowed
	so F(k) <= stress_allowed * A
	so F(k) <= (safety factor) * stress_yield * A 
These constant values are stored in "material.in"

In the calculations I've set up the input so that the total load applied is 1
This means that the forces we have so far (elements in Q) are fractions of the load
i.e. F(member k) = K * L, where K is the number and L is whatever load we apply
this means:
	F(k) = K * L <= (safety factor) * stress_yield * A
	so L <= (safety factor) * stress_yield * A / K
If we set L = [that mess] then we get the maximum load before the member yields

2) buckling
For buckling, only check the members under compression. Check that:
	F(k) <= pi^2 * E * I / l^2
As above we say F(k) = K * L, so the maximum load before the member buckles is:
	L = pi^2 * E * I / (l^2 * K)

The maximum load for the member we're looking at is the lower of the two limits.

The maximum load for the bridge is the lowest of all members' maximum loads. 
So we just go through the list we just calculated and find the lowest.


--- step 6 ---
now the results are just presented in a legible way. 
it is a horrible, horrible section of code
