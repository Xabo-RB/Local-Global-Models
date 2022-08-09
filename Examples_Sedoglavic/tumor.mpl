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


r1 := -(k3+k7)*x1+k4*x2:
r2 := k3*x1-(k4+a*k5+b*d*k5)*x2+k6*x3+k6*x4+k5*x2*x3+k5*x2*x4:
r3 :=  a*k5*x2-k6*x3-k5*x2*x4:
r4 := b*d*k5*x2-k6*x4-k5*x2*x4:
r5 := k7*x1:

VectorField := [

	r1,
	r2,
	r3,
	r4,
	r5
] :

# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

	x5
] :

#-----------------------------------------------------------------------------#
OutputsVariables:= [y1] 					:
Inputs 		:= [] 					:
Parameters 	:= [k3, k4, k5, k6, k7, a,b,d]					:
# The variables have to be ordered as the vectors field.
Variables 	:= [x1,x2,x3,x4,x5] 					:
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
