"""
    pin(allsteps; parent::Int = 0)

Add a pin 📌 for the given index
"""
function pin(allsteps; index::Int=0)
    index = _validindex(allsteps, index)
    allsteps.pinned[index] = true
    return allsteps
end

"""
    unpin(allsteps; index::Integer)

Remeve the pin 📌 with the given index
"""
function unpin(allsteps; index::Int=0)
    index = _validindex(allsteps, index)
    allsteps.pinned[index] = false
    return allsteps
end

function _find_all_pinned_states(allsteps)
    return findall(allsteps.pinned)
end
