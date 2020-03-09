struct TreeNode
    parent::Int
    children::Vector{Int}
    comment::String
end

struct Tree
    nodes::Vector{TreeNode}
end

Tree() = Tree([TreeNode(0, Vector{Int}(),"Initial Configuration")])

function addchild(tree::Tree, index::Int, comment::String)
    if index == 0
        index = length(tree.nodes)
    elseif index < 0
        throw(BoundsError(tree, index))
    end
    push!(tree.nodes, TreeNode(index, Vector{}(),comment))
    child = length(tree.nodes)
    push!(tree.nodes[index].children, child)
    child
end

children(tree, id) = tree.nodes[id].children
parent(tree,id) = tree.nodes[id].parent
comment(tree,id)=tree.nodes[id].comment
