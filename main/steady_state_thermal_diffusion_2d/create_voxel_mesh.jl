
include("node.jl")
include("element.jl")
#----------------------------------------------------------------
# モデルの長さ、分割数からボクセルメッシュを作成
# 2D version
#----------------------------------------------------------------
function create_voxel_mesh(length_x, length_y, division_x, division_y)

    # 節点の生成
    # 初期化
    num_node = (division_x + 1) * (division_y + 1)
    nodes = Vector{Node}(undef, num_node)
    
    # 節点間距離の計算
    delta_x::Float64 = length_x / division_x
    delta_y::Float64 = length_y / division_y
    
    # 分割数のループ
    for j = 1 : division_y+1
        # y座標
        coord_y::Float64 = (j - 1) * delta_y
        
        for i = 1 : division_x+1
            # x座標
            coord_x::Float64 = (i - 1) * delta_x
            
            # 節点番号
            id::Int64 = i + (j - 1) + division_x * (j - 1);

            # 情報の更新
            nodes[id] = Node(id, [coord_x, coord_y, 0.0], 0.0) 
        end
    end

    # コネクティビティの作成
    num_element = division_x * division_y
    connects = Matrix{Int64}(undef, num_element, 4)
    for j = 1 : division_y
        for i = 1 : division_x
            id = i + (j - 1) * division_x
            node_id1 = i + (j - 1) * (division_x + 1)
            node_id2 = node_id1 + 1
            node_id3 = node_id1 + division_x + 1 + 1
            node_id4 = node_id1 + division_x + 1
            connects[id, :] = [node_id1 node_id2 node_id3 node_id4]
        end
    end

    return num_node, num_element, nodes, connects
end
