function auto_allsteps(D::Tuple{ClosedInterval{<:Real}, ClosedInterval{<:Real}})
    auto_allsteps!(StepTree(), D)
end

function auto_allsteps(Ds::Vector{<:Tuple{ClosedInterval{<:Real}, ClosedInterval{<:Real}}})
    steptree = StepTree()
    for D in Ds
        auto_allsteps!(steptree, D)
    end
    return steptree
end

function auto_allsteps!(steptree::StepTree, D::Tuple{ClosedInterval{<:Real}, ClosedInterval{<:Real}})
    _, D₂ = D
    steptree = initial_state!(steptree, D)
    newton_onestep!(steptree, fixingmethod=:fix3points)
    newton_onestep!(steptree)
    refinement!(steptree, p₊=(0,1), k₊=(EmptyKnotVector(),KnotVector([mean(D₂)])))
    newton_onestep!(steptree)
    newton_onestep!(steptree)
    newton_onestep!(steptree)
    pin!(steptree)
    return steptree
end
