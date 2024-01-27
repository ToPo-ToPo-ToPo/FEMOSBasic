
include("node.jl")
include("element.jl")
using WriteVTK
#----------------------------------------------------------------
# Paraviewのvtkファイルに結果を書き出す関数
#----------------------------------------------------------------
function write_vtk_unstructured(nodes, elements, results, filename)
    
    # 節点情報の設定
    num_point = length(nodes) 
    points = Matrix{Float64}(undef, 3, num_point)
    for i in 1 : 3
        for node in nodes
            id = node.id
            points[i, id] = node.coordinate[i]
        end
    end

    # 要素タイプの設定
    num_cell = length(elements)
    cells = Vector{MeshCell}(undef, num_cell)
    for element in elements
        id = element.id
        node_ids = [element.nodes[1].id, element.nodes[2].id, element.nodes[3].id, element.nodes[4].id]
        cells[id] = MeshCell(VTKCellTypes.VTK_QUAD, node_ids)
    end

    # 結果の出力
    vtk_grid(filename, points, cells) do vtk
        vtk["temperature", VTKPointData()] = results
    end
    
    

end