
using SparseArrays
#----------------------------------------------------------------
# 材料モデル
#----------------------------------------------------------------
mutable struct Material
    conductivity::Float64           # 熱伝導率
end
#----------------------------------------------------------------
# 積分点構造体
#----------------------------------------------------------------
mutable struct EvaluatePoint
    id::Int64                       # 評価点の番号
    material::Material              # 材料モデル
end
#----------------------------------------------------------------
# 節点構造体
#----------------------------------------------------------------
mutable struct Node
    id::Int64                       # 節点番号
    coordinate::Vector{Float64}     # 節点座標
    field::Float64                  # 状態場 (温度)
end
#----------------------------------------------------------------
# 要素構造体
#----------------------------------------------------------------
mutable struct Element
    id::Int64                                # 要素番号
    nodes::Vector{Node}                      # 要素節点
    evaluate_points::Vector{EvaluatePoint}   # 積分点の情報

    #-------------------------------------------------------
    # 内部コンストラクタ
    #-------------------------------------------------------
    function Element(i::Int64, n::Vector{Node}, m::Material)
        # 要素番号の作成
        id = i

        # 要素内節点の作成
        nodes = n

        # 評価点の作成
        evaluate_points = [
            EvaluatePoint(1, deepcopy(m)),
            EvaluatePoint(2, deepcopy(m)),
            EvaluatePoint(3, deepcopy(m)),
            EvaluatePoint(4, deepcopy(m))
        ]

        return new(id, nodes, evaluate_points)
    end
end
#----------------------------------------------------------------
# 要素中心座標を計算する関数
#----------------------------------------------------------------
function compute_element_center(element::Element)

    # 初期化
    center_coord = zeros(Float64, 3)

    # 主な処理
    for node in element.nodes
        center_coord += node.coordinate
    end

    # 平均化した値を返す
    return center_coord / length(element.nodes)
end
#----------------------------------------------------------------
# ディレクレ境界条件の構造体
#----------------------------------------------------------------
mutable struct Dirichlet
    id::Int64                                 # 境界条件の番号
    value::Float64                            # 境界条件値
    nodes::Vector{Node}                       # 条件を与える節点
end
#----------------------------------------------------------------
# ノイマン境界条件の構造体
#----------------------------------------------------------------
mutable struct Neumman
    id::Int64                                 # 境界条件の番号
    value::Float64                            # 境界条件値
    nodes::Vector{Node}                       # 条件を与える節点
end
#----------------------------------------------------------------
# ソース項の構造体
#----------------------------------------------------------------
mutable struct Source
    id::Int64                                 # 境界条件の番号
    value::Float64                            # 境界条件値
    elements::Vector{Element}                 # 条件を与える要素
end
#----------------------------------------------------------------
# モデル
#----------------------------------------------------------------
mutable struct Model
    nodes::Vector{Node}                       # 全節点の情報
    elements::Vector{Element}                 # 全要素の情報
end
#----------------------------------------------------------------
# 有限要素法で2次元熱伝導方程式を解く関数
#----------------------------------------------------------------
function solve(model::Model, dirichlet_bcs::Vector{Dirichlet}, neumman_bcs::Vector{Neumman}, sources::Vector{Source})
    
    # 基本データの取得
    nodes = model.nodes
    elements = model.elements
    
    # 総節点数の取得
    num_nodes = length(nodes)

    # グローバル行列と全体ベクトルの初期化
    K  = sparse(zeros(Float64, num_nodes, num_nodes))
    Fq = zeros(Float64, num_nodes)
    Fb = zeros(Float64, num_nodes)

    # 全体熱流束ベクトルを作成
    for neumman_bc in neumman_bcs
        # 熱流束の値を取得
        value = neumman_bc.value
        
        # ベクトルに足し込む
        for node in neumman_bc.nodes
            Fq[node.id] += value
        end
    end

    # 全体熱伝導マトリクスを作成
    for element in elements
        # 要素行列作成
        #Ke = compute_Ke(element)

        # グローバル剛性行列と右辺ベクトルに組み込み
        #for i in 1:4
        #    for j in 1:4
        #        K_global[element.nodes[i], element.nodes[j]] += Ke[i, j]
        #    end
        #    F[element.nodes[i]] += local_F[i]
        #end
    end

    # ソース項ベクトルを作成
    for source in sources
        for element in source.elements
            # 要素ベクトルを計算
            #Fe = compute_Fbe(element)

            # アセンブリング
            #for i in 1:4
                #Fb[element.nodes[i]] += Fe[i]
            #end
        end
    end

    # 方程式の作成
    #lhs = K
    #rhs = Fq + Fb

    # 温度固定境界条件の適用
    #for (node_idx, temperature) in fixed_boundaries
    #    lhs[node_idx, :] .= 0.0
    #    lhs[node_idx, node_idx] = 1.0
    #    rhs[node_idx] = temperature
    #end

    # 疎行列を使って方程式を解く
    #U = lhs \ rhs
    U = Fq

    # 結果をノードに反映
    #for i in 1:num_nodes
    #    nodes[i].temperature = U[i]
    #end

    return U
end
#----------------------------------------------------------------
#----------------------------------------------------------------
function compute_element_matrices(nodes, element)
    num_nodes = length(element.nodes)
    conductivity = element.conductivity

    # ガウス積分の重みと積分点の座標
    gauss_weights = [1.0, 1.0]
    gauss_points = [
        (-1 / sqrt(3), -1 / sqrt(3)),
        (1 / sqrt(3), -1 / sqrt(3)),
        (1 / sqrt(3), 1 / sqrt(3)),
        (-1 / sqrt(3), 1 / sqrt(3))
    ]

    # 初期化
    local_K = zeros(Float64, num_nodes, num_nodes)
    local_F = zeros(Float64, num_nodes)

    for gauss_point in gauss_points
        xi, eta = gauss_point

        # 形状関数とヤコビアン行列の計算
        N = shape_functions(xi, eta)
        J = jacobian_matrix(nodes, element, xi, eta)

        # ガウス積分の重みとヤコビアンのスケーリング
        weight = gauss_weights[1] * gauss_weights[2] * det(J)

        # ローカル剛性行列の計算
        local_K += conductivity * (inv(J)' * inv(J)) * N * N' * weight

        # ローカル右辺ベクトルの計算（内部発熱項も含む）
        local_F += N * element.heat_source * weight
    end

    return local_K, local_F
end
#----------------------------------------------------------------
# 形状関数の定義
#----------------------------------------------------------------
function shape_functions(xi, eta)
    N = [0.25 * (1 - xi) * (1 - eta),
         0.25 * (1 + xi) * (1 - eta),
         0.25 * (1 + xi) * (1 + eta),
         0.25 * (1 - xi) * (1 + eta)]
    return N
end
#----------------------------------------------------------------
# ヤコビアン行列の計算
#----------------------------------------------------------------
function jacobian_matrix(nodes, element, xi, eta)
    num_nodes = length(element.nodes)
    dN_dxi, dN_deta = shape_function_derivatives(xi, eta)

    J = zeros(2, 2)

    for i in 1:num_nodes
        J += [nodes[element.nodes[i]].x, nodes[element.nodes[i]].y] * dN_dxi[i, :]
    end

    return J
end
#----------------------------------------------------------------
# 形状関数の導関数の計算
#----------------------------------------------------------------
function shape_function_derivatives(xi::Float64, eta::Float64)
    # 
    dN_dxi = [-0.25 * (1 - eta), 0.25 * (1 - eta), 0.25 * (1 + eta), -0.25 * (1 + eta)]
    dN_deta = [-0.25 * (1 - xi), -0.25 * (1 + xi), 0.25 * (1 + xi), 0.25 * (1 - xi)]

    return [dN_dxi dN_deta]
end

#----------------------------------------------------------------
# メイン関数
#----------------------------------------------------------------
function main()
    # 材料モデルの定義
    kappa = 1.0
    material = Material(kappa)

    # 節点の定義
    num_node = 4
    nodes = Vector{Node}(undef, num_node)
    nodes[1] = Node(1, [0.0, 0.0, 0.0], 0.0)
    nodes[2] = Node(2, [1.0, 0.0, 0.0], 0.0)
    nodes[3] = Node(3, [1.0, 1.0, 0.0], 0.0)
    nodes[4] = Node(4, [0.0, 1.0, 0.0], 0.0)

    # 要素の定義
    num_elem = 1

    # コネクティビティの作成
    connect = Matrix{Node}(undef, num_elem, 4)
    connect[1, 1] = nodes[1]
    connect[1, 2] = nodes[2]
    connect[1, 3] = nodes[3]
    connect[1, 4] = nodes[4]

    # 要素の定義
    elements = Vector{Element}(undef, num_elem)
    for e = 1 : num_elem
        elements[e] = Element(1, connect[e, :], material)
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
    thermal_diff_2d_model = Model(nodes, elements)

    # FEM解析の実行
    U = solve(thermal_diff_2d_model, dirichlet_bcs, neumman_bcs, sources)

    # 結果の整理

    # 結果のポスト
end
#----------------------------------------------------------------
# main関数の呼び出し
#----------------------------------------------------------------
main()
