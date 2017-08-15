using FastaIO
using PyPlot

export NodeData
abstract type NodeData end

export EmptyNodeData
type EmptyNodeData <: NodeData

end

export MyNodeData
type MyNodeData <: NodeData
  sequence::AbstractString

  MyNodeData(sequence::AbstractString) = new(sequence)
end

export TreeNode
type TreeNode
  parent::Nullable{TreeNode}
  children::Array{TreeNode,1}
  branchlength::Float64
  name::AbstractString
  data::NodeData
  nodeindex::Int
  seqindex::Int

  TreeNode(branchlength::Float64, name::AbstractString) = new(Nullable{TreeNode}(),TreeNode[],branchlength,name, EmptyNodeData(),0,0)
  TreeNode(branchlength::Float64, name::AbstractString, data::NodeData) = new(Nullable{TreeNode}(),TreeNode[],branchlength,name,data,0,0)
end

Base.start(node::TreeNode) = 1
Base.done(node::TreeNode,state) = length(node.children) == state-1
Base.next(node::TreeNode,state) = node.children[state], state+1
Base.length(node::TreeNode) = length(node.children)

export roottree
function roottree(root::TreeNode,index::Int=1)
  newroot = TreeNode(0.0, "")
  child = splice!(root.children,index)
  child.parent = Nullable{TreeNode}()
  addchild(newroot,child)
  br = child.branchlength/2.0
  child.branchlength = br
  root.branchlength = br
  addchild(newroot,root)
  return newroot
end

export isroot
function isroot(node::TreeNode)
  return isnull(node.parent)
end

export isleafnode
function isleafnode(node::TreeNode)
  return length(node.children) == 0
end

export addchild
function addchild(parent::TreeNode, child::TreeNode)
  if isnull(child.parent)
    push!(parent.children, child)
    child.parent = Nullable{TreeNode}(parent)
  else
    if insubtree(parent, child)
      throw(ArgumentError("Cannot move node to a subtree of itself!"))
    else
      children = get(child.parent).children
      index = findfirst(children, child)
      deleteat!(children, index)
      child.parent = Nullable{TreeNode}(parent)
      push!(parent.children, child)
    end
  end
end

export insubtree
function insubtree(node::TreeNode, subtree::TreeNode)
  if subtree == node
    return true
  else
    for child in subtree
      if insubtree(node, child)
        return true
      end
    end
  end
  return false
end


function getnewickhelper(node::TreeNode)
  if length(node.children) == 0
    return string(node.name,":", node.branchlength)
  else
    ret = join(AbstractString[getnewickhelper(child) for child in node], ",")
    return string("(",ret,")",node.name,":", node.branchlength)
  end
end

export getnewick
function getnewick(node::TreeNode)
  return string(getnewickhelper(node),";")
end

function treefromnewickhelper(newick::AbstractString)
  startindex = 1
  endindex = length(newick)
  lm = match(r"[\)][^\)]*$", newick)

  tag = newick
  childstrings = AbstractString[]
  if lm != nothing
    tag = newick[lm.offset+1:end]
    childstring = newick[startindex+1:lm.offset-1]

    child = ""
    a = 0
    for i=1:length(childstring)
      if childstring[i] == '('
        a += 1
        child = string(child,childstring[i])
      elseif childstring[i] == ')'
        a -= 1
        child = string(child,childstring[i])
      elseif childstring[i] == ',' && a == 0
        push!(childstrings, child)
        child = ""
      else
        child = string(child,childstring[i])
      end
    end
    if child != ""
      push!(childstrings, child)
    end
  end
  spl = split(tag,":")
  name = ""
  branchlength = 0.0
  if length(spl) > 0
    name = strip(spl[1])
  end
  if length(spl) > 1
    branchlength = parse(Float64,spl[2])
  end

  return childstrings,name,branchlength
end

export gettreefromnewick
"""
  gettreefromnewick(newick)

Returns a TreeNode object from `newick` string.
"""
function gettreefromnewick(newick::AbstractString)
  childstrings,name,branchlength = treefromnewickhelper(rstrip(newick,';'))
  node = TreeNode(branchlength,name, EmptyNodeData())
  for childstring in childstrings
    addchild(node, gettreefromnewick(childstring))
  end
  return node
end

export prettyprint
function prettyprintstring(node::TreeNode, spaces::Int=0)
  ret = string(repeat("----",spaces),"+++++ ", node.name, "(", node.branchlength,")", "\n")
  for child in node
      ret = string(ret, prettyprintstring(child,spaces+1))
  end
  return ret
end

export ExtraNodeData
type ExtraNodeData <: NodeData
  nodeindex::Int
  seqindex::Int
  ExtraNodeData(seqindex::Int) = new(0,seqindex)
end

export annotatetree
function annotatetree(node::TreeNode, seqnametoindex::Dict{AbstractString,Int}, nodelist::Array{TreeNode,1}=TreeNode[])
  push!(nodelist,node)
  node.data = ExtraNodeData(get(seqnametoindex,node.name,0))
  node.data.nodeindex = length(nodelist)
  node.nodeindex = length(nodelist)
  node.seqindex = get(seqnametoindex,node.name,0)
  for childnode in node
    annotatetree(childnode,seqnametoindex,nodelist)
  end
end

export getnodelist
function getnodelist(node::TreeNode, nodelist::Array{TreeNode,1}=TreeNode[])
  push!(nodelist,node)
  for childnode in node
    getnodelist(childnode,nodelist)
  end
  return nodelist
end

function getpattern(data::Array{Float64,3}, node::TreeNode, col::Int, pattern::Array{Int8,1}=Int8[])
  if isleafnode(node)
    for b in data[node.seqindex,col,:]
      push!(pattern, Int8(b))
    end
  else
    for childnode in node
      getpattern(data,childnode,col,pattern)
    end
  end

  return pattern
end

function getqmatrix(distmat::Array{Float64, 2})
    n = size(distmat)[1]
    newmat = zeros(size(distmat))
    for i in 1:n
        for j in i+1:n
            newmat[i,j] = ((n-2)*distmat[i,j] - sum(distmat[i,:])) - sum(distmat[j,:])
        end
    end
    return newmat
end

function diagonal_ignore_indmin(matrix)
    min_val = Inf
    ind_min_x = 0
    ind_min_y = 0
    for x in 1:(size(matrix)[1])
        for y in 1:(size(matrix)[2])
            if x != y && matrix[x,y] < min_val
                min_val = matrix[x,y]
                ind_min_x = x
                ind_min_y = y
            end
        end
    end
    return (ind_min_x,ind_min_y)
end

function treebuild(sequences::Array{String,1}; distmat = nothing, namearr::Array{String,1} = ["orig$i" for i in 1:length(sequences)], disallow_negative_length=true)
    if distmat == nothing
      distmat = dist_matrix(sequences, sequences)
    end
    nodes = [TreeNode(0.0, namearr[i], MyNodeData(sequences[i])) for i in 1:length(sequences)]
    for (i,node) in enumerate(nodes)
        node.seqindex = i
    end
    newnodenum = 1
    ind1, ind2 = 0,0
    
    if length(namearr) != length(sequences)
        error("Incorrect number of names")
    end
        
    for i in 1:(length(sequences)-2)
        n = size(distmat)[1]
        
        qmat = getqmatrix(Array{Float64,2}(distmat))

        ind1, ind2 = diagonal_ignore_indmin(qmat)
        
        if (ind1 == ind2)
            error("Something broke: Fix qMatrix to ignore diagonal zeros")
        end
            
        newnode = TreeNode(0.0, "inferred$newnodenum")
        newnodenum += 1
        
        nodes[ind1].branchlength = ((1/2) * distmat[ind1, ind2]) + (1/(2*(n-2)))*(sum(distmat[ind1,:]) - sum(distmat[ind2,:]))
        nodes[ind2].branchlength = distmat[ind1, ind2] - nodes[ind1].branchlength
       
    if disallow_negative_length 
      if (nodes[ind1].branchlength < 0)
        nodes[ind2].branchlength += abs(nodes[ind1].branchlength)
        nodes[ind1].branchlength = 0
      elseif (nodes[ind2].branchlength <0)
        nodes[ind1].branchlength += abs(nodes[ind2].branchlength)
        nodes[ind2].branchlength = 0
      end 
    end

        addchild(newnode, nodes[ind1])
        addchild(newnode, nodes[ind2])
        
        deleteat!(nodes, sort([ind1, ind2]))
        push!(nodes, newnode)
        
        remaininds = filter(x -> (x != ind1 && x != ind2), 1:n)
        
        #This connects the last two nodes
        if length(nodes) == 2
            nodes[1].branchlength = distmat[remaininds[1], ind1] - nodes[ind1].branchlength
            addchild(nodes[2], nodes[1])
            deleteat!(nodes, 1)
            break
        end
            
        newrow = []
        for j in remaininds
            push!(newrow, (distmat[j, ind1] + distmat[j,ind2] - distmat[ind1, ind2])/2)
        end
        
        distmat = hcat(distmat[remaininds, remaininds], newrow)
        push!(newrow, 0)
        distmat = vcat(distmat, newrow')
    end
    return nodes[1]
end
            
function treedepth(node::TreeNode)
    if isleafnode(node)
        return 1
    else
        return maximum([treedepth(node.children[i]) for i in 1:length(node.children)]) + 1
    end
end

function getnonleaflist(node::TreeNode, nonleaflist::Array{TreeNode,1}=TreeNode[])
    if !isleafnode(node)
        push!(nonleaflist, node)
    end
    for childnode in node
        getnonleaflist(childnode, nonleaflist)
    end
    return nonleaflist
end

function getleaflist(node::TreeNode, leaflist::Array{TreeNode,1}=TreeNode[])
    if isleafnode(node)
        push!(leaflist, node)
    end
    for childnode in node
        getleaflist(childnode, leaflist)
    end
    return leaflist
end

function getdistfromroot(node::TreeNode)
    if isnull(node.parent)
        return 0
    else
        return getdistfromroot(get(node.parent)) + node.branchlength
    end
end

function drawtree(root::TreeNode; xStart::Float64=0.0, yStart::Float64=0.0, xDist::Float64 = 2.0, yDist::Float64 = 1.0, scaled::Bool = false, extend::Bool = false, 
        names::Bool = false, reversed::Bool = false, xEnd::Float64 = 0.0, bubbles::Bool = false, 
        bubble_color_vector::Vector{String}=["#000000" for i in 1:length(getleaflist(root))],
        name_color_vector::Vector{String}=["#000000" for i in 1:length(getleaflist(root))],
        branch_color_vector::Vector{String}=["#000000" for i in 1:length(getleaflist(root))]
        )
    
    levels = treedepth(root)
    nodes = getleaflist(root)
    
    #WIP: Figure out a way to draw a tree with one node, or some other alternative
    if(length(nodes) == 1)
        error("Can't draw a tree with only one node")
    end
            
    yposdict = Dict{TreeNode,Float64}()
    yposarr = linspace(yDist, -yDist, length(nodes))
    xposarr = [getdistfromroot(node) for node in nodes]
    scale = scaled ? abs(maximum(abs.(xposarr)) /xDist) : 1
    reversemult = reversed? -1 : 1
    xposarr = (xposarr ./ scale) .* reversemult
    
    if xEnd != 0.0
        xStart = xEnd + (-1 * ((maximum(abs.(xposarr))) * reversemult))
    end
        
    const extenddist = 2
    
    name_size = maximum([length(node.name) for node in nodes]) + 1
    
    #Populate posdict, and draw leaves
    for (index, leaf) in enumerate(nodes)
        yposdict[leaf] = yposarr[index]
        
        if (!isnull(leaf.parent))
            plot([xStart + xposarr[index], xStart + reversemult * getdistfromroot(get(leaf.parent))/scale], [yposdict[leaf], yposdict[leaf]], color=branch_color_vector[leaf.seqindex])
        end
        #nodeextend = (min(abs(abs(xposarr[index]) - maximum(abs.(xposarr))), abs(extenddist/scale)))
        nodeextend = 0

        if (names)
            if !reversed
                annotate(string(leaf.name, ["-" for i in 1:(name_size - length(leaf.name))]...), [xStart + reversemult * (maximum(abs.(xposarr))), yposdict[leaf]], horizontalalignment="left", verticalalignment="center", family="monospace", color=name_color_vector[leaf.seqindex])
            else
                annotate(string(["-" for i in 1:(name_size - length(leaf.name))]... , leaf.name), [xStart + reversemult * (maximum(abs.(xposarr))), yposdict[leaf]], horizontalalignment="right", verticalalignment="center", family="monospace", color=name_color_vector[leaf.seqindex])
            end
            
        end
        if (bubbles)
            scatter(xStart + xposarr[index], yposdict[leaf], c=bubble_color_vector[leaf.seqindex], zorder=20)
        end
        if (extend)
            plot([xposarr[index] + xStart + reversemult * nodeextend, xStart + reversemult * maximum(abs.(xposarr))], [yposdict[leaf], yposdict[leaf]], color = "0.8", alpha=0.15)
        end
    end
    
    waitingNodes = []
    for level in 2:levels
        #Filter null, and already considered nodes
        nodes = [node.parent for node in nodes]
        nodes = filter(!isnull, nodes)
        nodes = [get(node) for node in nodes]
        
        append!(nodes, waitingNodes)
        
        waitingNodes = [node for node in nodes if treedepth(node) != level]
        nodes = [node for node in nodes if treedepth(node) == level]
        
        nodes = filter(x -> !haskey(yposdict, x), nodes)
        if length(nodes) == 0
            break
        end
        
        for node in nodes
            childPosArr = [yposdict[child] for child in node]
            yposdict[node] = mean(childPosArr)
            
            xPos = xStart + (getdistfromroot(node)/scale) * reversemult
            plot([xPos, xPos], [maximum(childPosArr), minimum(childPosArr)], "k-")
            
            if isnull(node.parent)
                xEnd = xStart
            else
                xEnd = xStart + (reversemult * (getdistfromroot(get(node.parent))/scale))
            end

            plot([xPos, xEnd], [yposdict[node], yposdict[node]], "k-")
        end
    end
    
    return (xStart, xStart + reversemult * (maximum(abs.(xposarr))), name_size)
end 

function gettotaldist(tree1::TreeNode, tree2::TreeNode, flow; ignorefreq::Bool = false, exp = 1)            
    yposarr1 = linspace(1, -1, length(getleaflist(tree1)))
    yposarr2 = linspace(1, -1, length(getleaflist(tree2)))
    
    if ignorefreq
        map(x-> x != 0? 1.0 : 0.0, flow)
    end
    
    ordered_flow = flow[getorder(tree1),getorder(tree2)]
    
    totaldist = 0
    for i in 1:size(flow)[1]
        for j in 1:size(flow)[2]
            if ordered_flow[i,j] != 0
                totaldist += ordered_flow[i,j] * (abs(yposarr1[i] - yposarr2[j]) ^ exp)
            end
        end
    end

    return totaldist
end 

#Takes exp parameter but doesn't actually use this... this is just to allow compatibility
function countcrossings(tree1::TreeNode, tree2::TreeNode, flow; ignorefreq::Bool = false, exp = 1)            
    yposarr1 = linspace(1, -1, length(getleaflist(tree1)))
    yposarr2 = linspace(1, -1, length(getleaflist(tree2)))
    
    if ignorefreq
        map(x-> x != 0? 1.0 : 0.0, flow)
    end
    
    ordered_flow = flow[getorder(tree1),getorder(tree2)]
    
    #Stores triple of: left-end pos of segment, right-end pos of segment, weight
    arrtriples = []
    
    totaldist = 0
    for i in 1:size(flow)[1]
        for j in 1:size(flow)[2]
            if ordered_flow[i,j] != 0
                push!(arrtriples, (yposarr1[i], yposarr2[j], ordered_flow[i,j]))
            end
        end
    end
        
    for i in 1:length(arrtriples)
        for j in i+1:length(arrtriples)
            if ((arrtriples[i][1] < arrtriples[j][1]) && (arrtriples[i][2] > arrtriples[j][2])) || 
                ((arrtriples[i][1] > arrtriples[j][1]) && (arrtriples[i][2] < arrtriples[j][2]))
                totaldist += arrtriples[i][3] * arrtriples[j][3]
            end
        end
    end

    return totaldist
end

function shuffletree(tree)
    for node in getnodelist(tree)
        if length(node.children) != 0
            shuffle!(node.children)
        end
    end
    return tree
end

function ladderize(tree)
    newtree = deepcopy(tree)
    for node in getnodelist(newtree)
        if length(node.children) != 0
            sort!(node.children, lt= (x,y)->length(getnodelist(x)) < length(getnodelist(y)))
        end
    end
    return newtree
end

function drawtreeswithflow(tree1, tree2, flow; figsize=(20,10), scaled=false, names = false, bubbles = false, 
        tree1_bubble_color_vector::Vector{String}=["#000000" for i in 1:length(getleaflist(tree1))],
        tree1_name_color_vector::Vector{String}=["#000000" for i in 1:length(getleaflist(tree1))],
        tree1_branch_color_vector::Vector{String}=["#000000" for i in 1:length(getleaflist(tree1))],
        tree2_bubble_color_vector::Vector{String}=["#000000" for i in 1:length(getleaflist(tree2))],
        tree2_name_color_vector::Vector{String}=["#000000" for i in 1:length(getleaflist(tree2))], 
        tree2_branch_color_vector::Vector{String}=["#000000" for i in 1:length(getleaflist(tree1))],
        thickness::Float64 = 2.0)
    figure(figsize=figsize)
    inter_tree_ratio = 0.5                                                            
    yDist = 1.0           
    min_width = 0.05
    char_width = 1/8
        
    max_right_branch = maximum(getdistfromroot.(getnodelist(tree2)))
    x1_start, x1_end, name_size_1 = drawtree(tree1; extend=true, names = names, reversed = false, yDist = yDist, bubble_color_vector=tree1_bubble_color_vector, name_color_vector=tree1_name_color_vector, branch_color_vector=tree1_branch_color_vector, scaled=scaled, bubbles=bubbles)
    
    trees_size = (x1_end + max_right_branch)
    inter_tree_dist = inter_tree_ratio * trees_size
    total_size = trees_size + inter_tree_dist
    
    #points per inch  = points/inch = total_size/figsize[1]
    #word_size = name_size_1 * width of each character * points per inch 
    word_size_1 = (name_size_1 * char_width * ( total_size / figsize[1]))
    
    x2_start, x2_end, name_size_2 = drawtree(tree2; extend=true, names = names, reversed = true, xEnd = x1_end + inter_tree_dist, yDist = yDist, bubble_color_vector=tree2_bubble_color_vector, branch_color_vector=tree2_branch_color_vector, name_color_vector=tree2_name_color_vector, scaled=scaled, bubbles=bubbles)
                                                 
    word_size_2 = (name_size_2 * char_width * ( total_size / figsize[1]))
    
    y1arr = linspace(yDist, -yDist, size(flow)[1])                                  
    y2arr = linspace(yDist, -yDist, size(flow)[2])  
    
    order1 = getorder(tree1)
    order2 = getorder(tree2)
    
    ordered_flow = flow[order1,order2]
                                                                                    
    for i in 1:size(flow)[1]                                                        
        for j in 1:size(flow)[2]                                                    
            if ordered_flow[i,j] != 0                                                       
                plot([x1_end + word_size_1, x2_end - word_size_2], [y1arr[i], y2arr[j]], linewidth=thickness*(min_width + ordered_flow[i,j]), "g-")
            end                                                                     
        end                         
    end                             
    return x1_start,x2_start
end 

function getorder(tree::TreeNode)
    return [node.seqindex for node in getleaflist(tree)]
end

function randopttrees(tree1::TreeNode, tree2::TreeNode, flow; tries=100, lossfunc=gettotaldist, ladderize1::Bool = false, ladderize2::Bool = false, exp=1, print_scores::Bool = false)    
    warn("randopttrees destructively reorders trees")
    if ladderize1
        ladderize(tree1)
    end
    if ladderize2
        ladderize(tree2)
    end
    
    min1 = deepcopy(tree1)
    min2 = deepcopy(tree2)
        
    mindist = lossfunc(tree1, tree2, flow, exp=exp, ignorefreq=true)
    
    if print_scores
        println("Original loss: ", mindist)
    end
    minflow = flow
    
    tried = 0
    while tried < tries
        copy1 = deepcopy(tree1)
        copy2 = deepcopy(tree2)    
        
        #First coin flip to choose which tree to alter
        #Second coin flip to choose to either shuffle a subtree, or just swap its children
        if (ladderize2 || rand()) > 0.5 && !ladderize1 
            shuffletree(tree1)
        else
            shuffletree(tree2)
        end
       
        orderarr1 = getorder(tree1)
        orderarr2 = getorder(tree2)
                        
        if (dist = lossfunc(tree1, tree2, flow, exp=exp, ignorefreq=true)) < mindist
            min1 = deepcopy(tree1)
            min2 = deepcopy(tree2)
            mindist = dist
            tried += 1
        else 
            tried += 1
            tree1 = copy1
            tree2 = copy2
        end    
    end
    if print_scores
        println("Optimized loss: ", mindist)
    end
    return min1, min2
end

function greedyopttrees(tree1::TreeNode, tree2::TreeNode, flow;fails_till_stop=100, lossfunc=gettotaldist, ladderize1::Bool = false, ladderize2::Bool = false, exp=1, print_scores::Bool = false)
    warn("greedyopttrees destructively reorders trees")
    if ladderize1
        ladderize(tree1)
    end
    if ladderize2
        ladderize(tree2)
    end
    
    min1 = deepcopy(tree1)
    min2 = deepcopy(tree2)
        
    mindist = lossfunc(tree1, tree2, flow, exp=exp, ignorefreq=true)
    
    if print_scores
        println("Original loss: ", mindist)
    end
    minflow = flow
    
    successes = 0
    fails = 0
    while fails < fails_till_stop
        copy1 = deepcopy(tree1)
        copy2 = deepcopy(tree2)    
        
        #First coin flip to choose which tree to alter
        #Second coin flip to choose to either shuffle a subtree, or just swap its children
        if (ladderize2 || rand()) > 0.5 && !ladderize1 
            if rand() > 0.5
                shuffletree(sample(getnonleaflist(tree1)))
            else
                reverse!(sample(getnonleaflist(tree1)).children)
            end
        else
            if rand() > 0.5
                shuffletree(sample(getnonleaflist(tree2)))
            else
                reverse!(sample(getnonleaflist(tree2)).children)
            end
        end
       
        orderarr1 = getorder(tree1)
        orderarr2 = getorder(tree2)
                        
        if (dist = lossfunc(tree1, tree2, flow, exp=exp, ignorefreq=true)) < mindist
            min1 = deepcopy(tree1)
            min2 = deepcopy(tree2)
            mindist = dist
            fails = 0
            successes +=1 
        else 
            fails += 1
            tree1 = copy1
            tree2 = copy2
        end    
    end
    if print_scores
        println("Optimized loss: ", mindist)
        println("Changes made: ", successes)
    end
    return min1, min2
end

#Wrapper function that will contain calls to various optimization functions
function optimize_trees(tree1::TreeNode, tree2::TreeNode, flow; ladderize1::Bool=false, ladderize2::Bool=false, print_scores::Bool=false)
    warn("Ignore subsequent warnings: optimize_trees prevents destructive reordering")
    #Take a copy to prevent modifying original
    tree1 = deepcopy(tree1)
    tree2 = deepcopy(tree2)
    
    #Minimize sum of squares of the vertical distance
    tree1, tree2 = randopttrees((tree1), (tree2), flow, tries=1000, lossfunc=gettotaldist, exp=2, ladderize1=ladderize1, ladderize2 = ladderize2, print_scores=print_scores);
    tree1, tree2 = greedyopttrees(tree1, tree2, flow, fails_till_stop=5000, lossfunc=gettotaldist, exp=2, ladderize1=ladderize1, ladderize2 = ladderize2, print_scores=print_scores);
    
    #Minimize the total crossings of trees
    tree1, tree2 = randopttrees(tree1, tree2, flow, tries=1000, lossfunc=countcrossings, ladderize1=ladderize1, ladderize2 = ladderize2, print_scores=print_scores);
    tree1, tree2 = greedyopttrees(tree1, tree2, flow, fails_till_stop=5000, lossfunc=countcrossings, ladderize1=ladderize1, ladderize2 = ladderize2, print_scores=print_scores);
    
    return tree1, tree2
end 
