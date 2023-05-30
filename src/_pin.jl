"""
    pin(steptree, parent::Int = 0)

Add a pin 📌 for the given index
"""
function pin!(steptree, index::Int = 0)
    index = _validindex(steptree, index)
    steptree.pinned[index] = true
    return steptree
end

"""
    unpin(steptree, index::Integer)

Remeve the pin 📌 with the given index
"""
function unpin!(steptree, index::Int = 0)
    index = _validindex(steptree, index)
    steptree.pinned[index] = false
    return steptree
end
