module POV_Ray

using IntervalSets

export MESH2Export

function MESH2Export(p,D;name="MESH2Export.inc",mesh=(10,10),smooth=true,texture=nothing)
    POVvec(x)="<"*repr(x)[2:end-1]*"> "
    if(texture!=Nothing)
        POVtxt(x)="texture{pigment{"*texture[1]*POVvec(x)*"}}"
        POVvec2(x)=POVvec(x)*repr(x)[2:end-1]*" "
    end

    t=(collect(range(minimum(D[1]),stop=maximum(D[1]),length=mesh[1]+1)),collect(range(minimum(D[2]),stop=maximum(D[2]),length=mesh[2]+1)))
    tc=(collect(range(minimum(D[1])+(maximum(D[1])-minimum(D[1]))/(2*mesh[1]),stop=maximum(D[1])-(maximum(D[1])-minimum(D[1]))/(2*mesh[1]),length=mesh[1])),
        collect(range(minimum(D[2])+(maximum(D[2])-minimum(D[2]))/(2*mesh[2]),stop=maximum(D[2])-(maximum(D[2])-minimum(D[2]))/(2*mesh[2]),length=mesh[2])))

    P=[p([t1,t2])[i] for t1 in t[1], t2 in t[2], i in 1:3]
    Pc=[p([t1,t2])[i] for t1 in tc[1], t2 in tc[2], i in 1:3]
    if(smooth)
        ṗ(t)=ForwardDiff.jacobian(p,t)
        e3(t)=cross(ṗ(t)[1:3,1],ṗ(t)[1:3,2])#LinearAlgebraは不要??
        E=[e3([t1,t2])[i] for t1 in t[1], t2 in t[2], i in 1:3]
        Ec=[e3([t1,t2])[i] for t1 in tc[1], t2 in tc[2], i in 1:3]
    end
    if(texture!=nothing)
        colo=texture[2]
        T=[colo([t1,t2])[i] for t1 in t[1], t2 in t[2], i in 1:3]
        Tc=[colo([t1,t2])[i] for t1 in tc[1], t2 in tc[2], i in 1:3]
    end

    Fc=[i₁-1+mesh[1]*(i₂-1) for i₁ in 1:mesh[1], i₂ in 1:mesh[2]]
    F=[mesh[1]*mesh[2]+i₁-1+(mesh[1]+1)*(i₂-1) for i₁ in 1:(mesh[1]+1), i₂ in 1:(mesh[2]+1)]
    Fi₁=hcat(reshape(Fc,mesh[1]*mesh[2]),reshape(F[1:mesh[1],1:mesh[2]],mesh[1]*mesh[2]),reshape(F[2:(mesh[1]+1),1:mesh[2]],mesh[1]*mesh[2]))
    Fi₂=hcat(reshape(Fc,mesh[1]*mesh[2]),reshape(F[2:(mesh[1]+1),1:mesh[2]],mesh[1]*mesh[2]),reshape(F[2:(mesh[1]+1),2:(mesh[2]+1)],mesh[1]*mesh[2]))
    Fi₃=hcat(reshape(Fc,mesh[1]*mesh[2]),reshape(F[2:(mesh[1]+1),2:(mesh[2]+1)],mesh[1]*mesh[2]),reshape(F[1:mesh[1],2:(mesh[2]+1)],mesh[1]*mesh[2]))
    Fi₄=hcat(reshape(Fc,mesh[1]*mesh[2]),reshape(F[1:mesh[1],2:(mesh[2]+1)],mesh[1]*mesh[2]),reshape(F[1:mesh[1],1:mesh[2]],mesh[1]*mesh[2]))

    MESH2=
    "mesh2{\n"*
    "vertex_vectors{\n"*
    string((mesh[1]+1)*(mesh[2]+1)+mesh[1]*mesh[2])*"\n"*
    *(map(POVvec,reshape([Pc[i₁,i₂,:] for i₁ in 1:mesh[1], i₂ in 1:mesh[2]],mesh[1]*mesh[2]))...)*"\n"*
    *(map(POVvec,reshape([P[i₁,i₂,:] for i₁ in 1:(mesh[1]+1), i₂ in 1:(mesh[2]+1)],(mesh[1]+1)*(mesh[2]+1)))...)*"\n"*
    "}\n"
    if(smooth)
        MESH2=MESH2*
        "normal_vectors{\n"*
        string((mesh[1]+1)*(mesh[2]+1)+mesh[1]*mesh[2])*"\n"*
        *(map(POVvec,reshape([Ec[i₁,i₂,:] for i₁ in 1:mesh[1], i₂ in 1:mesh[2]],mesh[1]*mesh[2]))...)*"\n"*
        *(map(POVvec,reshape([E[i₁,i₂,:] for i₁ in 1:(mesh[1]+1), i₂ in 1:(mesh[2]+1)],(mesh[1]+1)*(mesh[2]+1)))...)*"\n"*
        "}\n"
    end
    if(texture!=nothing)
        MESH2=MESH2*
        "texture_list{\n"*
        string((mesh[1]+1)*(mesh[2]+1)+mesh[1]*mesh[2])*"\n"*
        *(map(POVtxt,reshape([Tc[i₁,i₂,:] for i₁ in 1:mesh[1], i₂ in 1:mesh[2]],mesh[1]*mesh[2]))...)*"\n"*
        *(map(POVtxt,reshape([T[i₁,i₂,:] for i₁ in 1:(mesh[1]+1), i₂ in 1:(mesh[2]+1)],(mesh[1]+1)*(mesh[2]+1)))...)*"\n"*
        "}\n"
        MESH2=MESH2*
        "face_indices{\n"*
        string(4*mesh[1]*mesh[2])*"\n"*
        *(map(POVvec2,[Fi₁[i,:] for i in 1:(mesh[1]*mesh[2])])...)*"\n"*
        *(map(POVvec2,[Fi₂[i,:] for i in 1:(mesh[1]*mesh[2])])...)*"\n"*
        *(map(POVvec2,[Fi₃[i,:] for i in 1:(mesh[1]*mesh[2])])...)*"\n"*
        *(map(POVvec2,[Fi₄[i,:] for i in 1:(mesh[1]*mesh[2])])...)*"\n"*
        "}\n"*
        "}"
    else
        MESH2=MESH2*
        "face_indices{\n"*
        string(4*mesh[1]*mesh[2])*"\n"*
        *(map(POVvec,[Fi₁[i,:] for i in 1:(mesh[1]*mesh[2])])...)*"\n"*
        *(map(POVvec,[Fi₂[i,:] for i in 1:(mesh[1]*mesh[2])])...)*"\n"*
        *(map(POVvec,[Fi₃[i,:] for i in 1:(mesh[1]*mesh[2])])...)*"\n"*
        *(map(POVvec,[Fi₄[i,:] for i in 1:(mesh[1]*mesh[2])])...)*"\n"*
        "}\n"*
        "}"
    end

    fp = open(name,"w")
    write(fp, MESH2)
    close(fp)
    return true
end

end
