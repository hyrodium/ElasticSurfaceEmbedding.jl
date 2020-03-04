struct TreeNode
    parent::Int
    children::Vector{Int}
    comment::String
end

struct Tree
    nodes::Vector{TreeNode}
end

Tree() = Tree([TreeNode(0, Vector{Int}(),"Initial Configuration")])

function addchild(tree::Tree, id::Int, comment::String)
    1 <= id <= length(tree.nodes) || throw(BoundsError(tree, id))
    push!(tree.nodes, TreeNode(id, Vector{}(),comment))
    child = length(tree.nodes)
    push!(tree.nodes[id].children, child)
    child
end

children(tree, id) = tree.nodes[id].children
parent(tree,id) = tree.nodes[id].parent
comment(tree,id)=tree.nodes[id].comment

function shownode(tree,id,depth)
    txt="  "^depth*string(id)*": "*comment(tree,id)*"\n"
    for node ∈ children(tree,id)
        txt=txt*shownode(tree,node,depth+1)
    end
    return txt
end

function showtree(tree)
    shownode(tree,1,0)
end
