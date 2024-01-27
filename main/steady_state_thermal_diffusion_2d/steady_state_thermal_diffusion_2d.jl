
include("model.jl")
include("create_voxel_mesh.jl")
#----------------------------------------------------------------
# メイン関数
#----------------------------------------------------------------
function main()
    # 材料モデルの定義
    kappa = 1.0
    material = Material(kappa)

    # 解析モデルの設定
    length_x = 1.0
    length_y = 1.0
    division_x = 100
    division_y = 100

    # 解析モデルの作成
    num_node, num_element, nodes, connects = create_voxel_mesh(length_x, length_y, division_x, division_y)

    # 要素の定義
    elements::Vector{Element} = Vector{Element}(undef, num_element)
    for e = 1 : length(elements)
        nodeset::Vector{Node} = [nodes[connects[e, 1]], nodes[connects[e, 2]], nodes[connects[e, 3]], nodes[connects[e, 4]]]
        elements[e] = Element(e, nodeset, material)
    end

    # 節点の定義
    #num_node::Int64 = 4
    #nodes::Vector{Node} = Vector{Node}(undef, num_node)
    #nodes[1] = Node(1, [0.0, 0.0, 0.0], 0.0)
    #nodes[2] = Node(2, [1.0, 0.0, 0.0], 0.0)
    #nodes[3] = Node(3, [1.0, 1.0, 0.0], 0.0)
    #nodes[4] = Node(4, [0.0, 1.0, 0.0], 0.0)

    # 要素の定義
    #num_elem::Int64 = 1

    # コネクティビティの作成
    #connect = Matrix{Int64}(undef, num_elem, 4)
    #connect[1, :] = [1 2 3 4]

    # 要素の定義
    #elements::Vector{Element} = Vector{Element}(undef, length(connects))
    #for e = 1 : num_elem
    #    nodeset::Vector{Node} = [nodes[connect[e, 1]], nodes[connect[e, 2]], nodes[connect[e, 3]], nodes[connect[e, 4]]]
    #    elements[e] = Element(e, nodeset, material)
    #end

    # ディレクレ境界条件を設定する節点ベクトルを作成
    range_min = [0.0, 0.0, 0.0]
    range_max = [0.0, 1.0, 0.0]
    dirich_node = Vector{Node}()
    for node in nodes
        if range_min[1] - 1.0e-05 < node.coordinate[1] < range_max[1] + 1.0e-05 && 
           range_min[2] - 1.0e-05 < node.coordinate[2] < range_max[2] + 1.0e-05 &&
           range_min[3] - 1.0e-05 < node.coordinate[3] < range_max[3] + 1.0e-05
            
            # 該当する節点を追加
            push!(dirich_node, node)
        end
    end

    # ディレクレ境界の作成
    dirichlet_bcs = Vector{Dirichlet}()
    push!(dirichlet_bcs, Dirichlet(1, 0.0, dirich_node))

    # ノイマン境界を設定する節点ベクトルを作成
    range_min = [1.0, 0.0, 0.0]
    range_max = [1.0, 1.0, 0.0]
    neum_node = Vector{Node}()
    for node in nodes
        if range_min[1] - 1.0e-05 < node.coordinate[1] < range_max[1] + 1.0e-05 && 
           range_min[2] - 1.0e-05 < node.coordinate[2] < range_max[2] + 1.0e-05 &&
           range_min[3] - 1.0e-05 < node.coordinate[3] < range_max[3] + 1.0e-05
            
            # 該当する節点を追加
            push!(neum_node, node)
        end
    end

    # ノイマン境界の作成
    neumman_bcs = Vector{Neumman}()
    push!(neumman_bcs, Neumman(1, 1.0, neum_node))

    # ソース項を与える要素ベクトルを作成
    range_min = [0.0, 0.0, 0.0]
    range_max = [1.0, 1.0, 0.0]
    source_element = Vector{Element}()
    for element in elements
        # 要素中心を取得
        center_coord = compute_element_center(element)
        
        # 含まれるかどうか確認
        if range_min[1] - 1.0e-05 < center_coord[1] < range_max[1] + 1.0e-05 && 
           range_min[2] - 1.0e-05 < center_coord[2] < range_max[2] + 1.0e-05 &&
           range_min[3] - 1.0e-05 < center_coord[3] < range_max[3] + 1.0e-05
            
            # 該当する節点を追加
            push!(source_element, element)
        end
    end

    # ソース項の作成
    sources = Vector{Source}()
    push!(sources, Source(1, 0.0, source_element))

    # モデルの作成
    thermal_diff_2d_model::Model = Model(nodes, elements)

    # FEM解析の実行
    T = solve(thermal_diff_2d_model, dirichlet_bcs, neumman_bcs, sources)

    # 結果をノードに反映
    for node in thermal_diff_2d_model.nodes
        node.T = T[node.id]
    end

    # 結果のポスト
    for node in thermal_diff_2d_model.nodes
        println("node ID: $(node.id), (x, y) = ($(node.coordinate[1]), $(node.coordinate[2])), temperature = $(node.T)")
    end

end
#----------------------------------------------------------------
# main関数の呼び出し
#----------------------------------------------------------------
main()
