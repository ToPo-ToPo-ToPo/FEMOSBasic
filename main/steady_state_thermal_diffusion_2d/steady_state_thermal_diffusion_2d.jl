
include("model.jl")
include("create_voxel_mesh.jl")
include("output_vtu.jl")
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
    division_x = 200
    division_y = 200

    # 解析モデルの作成
    num_node, num_element, nodes, connects = create_voxel_mesh(length_x, length_y, division_x, division_y)

    # 要素の定義
    elements::Vector{Element} = Vector{Element}(undef, num_element)
    for e = 1 : length(elements)
        nodeset::Vector{Node} = [nodes[connects[e, 1]], nodes[connects[e, 2]], nodes[connects[e, 3]], nodes[connects[e, 4]]]
        elements[e] = Element(e, nodeset, material)
    end

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
    range_min = [1.0, 0.4, 0.0]
    range_max = [1.0, 0.6, 0.0]
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

    # vtkファイルに結果を書き出し
    write_vtk_unstructured(nodes, elements, T, "output/vtu_file/steady_state_thermal_diffusion")

end
#----------------------------------------------------------------
# main関数の呼び出し
#----------------------------------------------------------------
@time main()
