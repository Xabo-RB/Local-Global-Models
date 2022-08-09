#-----------------------------------------------------------------------------#
#								     	System 
# Description 	: 
# Result	:
#
#-----------------------------------------------------------------------------#
infolevel[observabilityTest] := 1 :
infolevel[observabilitySymmetries] := 1 :
#-----------------------------------------------------------------------------#
# Bibliography : Milena Anguelova, private communication 2006
# This examples is observable while previous code say opposite. The bug is know
# fixed. 
#-----------------------------------------------------------------------------#
# We assume that diff(Variables[i],t) = VectorsField[i]


r1 := -bi*x1*x2 - x1*mu - x1*bw*x3 + x4*a + mu:
r2 := bi*x1*x2 - gam*x2 + x1*bw*x3 - mu*x2:
r3 := xi*x2 - xi*x3:
r4 := gam*x2 - x4*mu - x4*a:

VectorField := [

	r1,
	r2,
	r3,
	r4
] :

# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

	k*x2,x1 + x4 + x2
] :

#-----------------------------------------------------------------------------#
OutputsVariables:= [y1,y2] 					:
Inputs 		:= [] 					:
Parameters 	:= [bi, gam, mu, bw, k, xi, a]					:
# The variables have to be ordered as the vectors field.
Variables 	:= [x1,x2,x3,x4] 					:
#-----------------------------------------------------------------------------#
# CAUTION read section II of the file INSTALL
libname := cat(currentdir(),"/release"),libname :
readlib(observabilityTest) :
NonObservable := observabilityTest(	VectorField	,
					Variables	,
					OutputSystem	,
					Parameters	,
					Inputs			) :
print(%) :					
GroupInfGen := observabilitySymmetries(	VectorField	,
					Variables	,
					OutputSystem	,
					Parameters	,
					Inputs		,
					NonObservable		) :
print(%) :
#-----------------------------------------------------------------------------#
quit :
