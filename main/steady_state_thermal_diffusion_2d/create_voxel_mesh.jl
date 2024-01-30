
include("node.jl")
include("element.jl")
include("line.jl")
#----------------------------------------------------------------
# モデルの長さ、分割数からボクセルメッシュを作成
# 2D version
#----------------------------------------------------------------
function create_voxel_mesh(length_x, length_y, division_x, division_y, thickness, material)

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

     # 要素の定義
     num_element = division_x * division_y
     elements = Vector{Element}(undef, num_element)
     for j = 1 : division_y
         for i = 1 : division_x
             id = i + (j - 1) * division_x
             n1 = i + (j - 1) * (division_x + 1)
             n2 = n1 + 1
             n3 = n1 + division_x + 1 + 1
             n4 = n1 + division_x + 1
             elements[id] = Element(id, [nodes[n1], nodes[n2], nodes[n3], nodes[n4]], thickness, material)
         end
     end
 

    # 線要素の定義
    num_line = division_x * (division_y + 1) + division_y * (division_x + 1)
    lines = Vector{Line}(undef, num_line)
    # X軸方向の線要素作成
    for j in 1 : division_y + 1
        for i in 1 : division_x
            id = (j - 1) * division_x + i
            n1 = i + (j - 1) * (division_x + 1)
            n2 = n1 + 1
            lines[id] = Line(id, [nodes[n1], nodes[n2]], thickness)
        end
    end

    # Y軸方向の線要素作成
    for j in 1 : division_y
        for i in 1 : division_x + 1
            id = division_x * (division_y + 1) + (i - 1) * division_y + j
            n1 = i + (j - 1) * (division_x + 1)
            n2 = n1 + division_x + 1
            lines[id] = Line(id, [nodes[n1], nodes[n2]], thickness)
        end
    end

    return nodes, elements, lines
end
