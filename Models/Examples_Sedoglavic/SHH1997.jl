using StructuralIdentifiability

ode = @ODEmodel(
    X'(t) = kcatX*X(t)*RVV/(kmX+X(t)),
    Xa'(t) = kcatX*X(t)*RVV/(kmX+X(t)) - kiXa*Xa(t) - kPT*Va(t)*Xa(t)*PL(t) + kPL*PT(t),
    V'(t) = -kcatV*V(t)*IIa(t)/(kmV+V(t)),
    Va'(t) = kcatV*V(t)*IIa(t)/(kmV+V(t)) - kPT*Va(t)*Xa(t)*PL(t) + kPL*PT(t),
    PL'(t) = - kPT*Va(t)*Xa(t)*PL(t) + kPL*PT(t),
    PT'(t) =  kPT*Va(t)*Xa(t)*PL(t) - kPL*PT(t),
    II'(t) = - kcatII*II(t)*PT(t)/(kmII+II(t)) - kcat2*II(t)*Xa(t)/(km2+II(t)),
    IIa'(t) =  kcatII*II(t)*PT(t)/(kmII+II(t)) + kcat2*II(t)*Xa(t)/(km2+II(t)) - kiIIaa2M*IIa(t) - kiIIaATIII*IIa(t),
    IIaa2M'(t) = kiIIaATIII*IIa(t),
    y1(t) = IIa(t) + 556*IIaa2M(t)/1000
)

@time println(assess_identifiability(ode))

