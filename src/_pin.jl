"""
    pin(; parent::Int = 0)

Add a pin 📌 for the given index
"""
function pin(allsteps; index::Int=0)
    index = _realparent(allsteps, index)
    allsteps.steps[index][1].pinned = true
    return allsteps
end

"""
    unpin(index::Integer)

Remeve the pin 📌 with the given index
"""
function unpin(allsteps; index::Int=0)
    index = _realparent(allsteps, index)
    allsteps.steps[index][1].pinned = false
    return allsteps
end

function _find_all_pinned_states(allsteps)
    return findall(step -> step[1].pinned, allsteps.steps)
end
