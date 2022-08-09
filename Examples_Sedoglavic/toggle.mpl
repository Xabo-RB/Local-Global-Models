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


r1 := k01 + k1/(1+(x2/(1+(u1/tatc)^natc))^ntetr) - x1:
r2 := k02 + k2/(1+(x1/(1+(u2/tiptg)^niptg))^nlaci) - x2:

VectorField := [

	r1,
	r2
] :

# We assume that OutputsVariables[i] = OutputSystem[i].
OutputSystem := [

	x1,x2
] :

#-----------------------------------------------------------------------------#
OutputsVariables:= [y1,y2] 					:
Inputs 		:= [u1,u2] 					:
Parameters 	:= [k01,k1,tatc,natc,ntetr,k02,k2,tiptg,niptg,nlaci]					:
# The variables have to be ordered as the vectors field.
Variables 	:= [x1,x2] 					:
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
