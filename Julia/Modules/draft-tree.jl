tree=Dict("1" => Dict("parent"=>"0","comment"=>"hogehoge1"),
          "2" => Dict("parent"=>"1","comment"=>"hogehoge2"),
          "3" => Dict("parent"=>"2","comment"=>"hogehoge3"),
          "4" => Dict("parent"=>"0","comment"=>"hogehoge4"),
          "5" => Dict("parent"=>"0","comment"=>"hogehoge5"),
          "6" => Dict("parent"=>"4","comment"=>"hogehoge6"),
          "7" => Dict("parent"=>"4","comment"=>"hogehoge7"),
          "8" => Dict("parent"=>"5","comment"=>"hogehoge8"),
          "9" => Dict("parent"=>"8","comment"=>"hogehoge9"),
          "10" => Dict("parent"=>"1","comment"=>"hogehoge10"),
          "11" => Dict("parent"=>"0","comment"=>"hogehoge11"),
          "12" => Dict("parent"=>"10","comment"=>"hoge12"))

function NodeSeries(tree,node)
    Nodes=[node]
    while Nodes[end] ≠ "0"
        push!(Nodes, tree[Nodes[end]]["parent"])
    end
    return Nodes
end

function TreeString(tree)
    serieses=Array{Int,1}[]
    for key in keys(tree)
        push!(serieses,(s->parse(Int,s)).(reverse(NodeSeries(tree,key))))
    end
    sort!(serieses)
    lowstrings=String[]
    n = length(serieses)
    for i in 1:n
        l=length(serieses[i])
        key=string(serieses[i][end])
        comment=tree[key]["comment"]
        if l == 2
            lowstring=key*": "*comment
            push!(lowstrings,lowstring)
        elseif l ≥ 3
            lowstring="  "^(l-3)*"└─"*key*": "*comment
            push!(lowstrings,lowstring)
            for j in 1:(i-1)
                chars=collect(lowstrings[end-j])
                if chars[2(l-3)+1]==' '
                    lowstrings[end-j]=join(chars[1:2(l-3)])*"│"*join(chars[2(l-3)+2:end])
                elseif chars[2(l-3)+1]=='└'
                    lowstrings[end-j]=join(chars[1:2(l-3)])*"├"*join(chars[2(l-3)+2:end])
                    break
                else
                    break
                end
            end
        end
    end
    outsting=""
    for s in lowstrings
        outsting=outsting*s*"\n"
    end
    return outsting
end

print(TreeString(tree))



println()

using JSON
open(ElasticSurfaceEmbedding.DIR*"/"*ElasticSurfaceEmbedding.NAME*".json","r") do f
    global dict=JSON.parse(f)
end

dict["2"]

print(TreeString(dict))

NodeSeries(dict,"7")

keys(dict)

delete!(dict,"Expr")
