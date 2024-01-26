
include("node.jl")
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
    coordinate::Vector{Float64}     # 積分点座標
    weight::Float64                 # 積分点の重み
    material::Material              # 材料モデル
end
#----------------------------------------------------------------
# 要素の構造体
#----------------------------------------------------------------
mutable struct Element
    id::Int64                                # 要素番号
    nodes::Vector{Node}                      # 要素節点
    evaluate_points::Vector{EvaluatePoint}   # 積分点の情報

    #-------------------------------------------------------
    # 内部コンストラクタ
    #-------------------------------------------------------
    function Element(id::Int64, nodes::Vector{Node}, material::Material)
        
        # 評価点の作成
        evaluate_points = [
            EvaluatePoint(1, [-1.0 / sqrt(3.0), -1.0 / sqrt(3.0)], 1.0, deepcopy(material)),
            EvaluatePoint(2, [ 1.0 / sqrt(3.0), -1.0 / sqrt(3.0)], 1.0, deepcopy(material)),
            EvaluatePoint(3, [ 1.0 / sqrt(3.0),  1.0 / sqrt(3.0)], 1.0, deepcopy(material)),
            EvaluatePoint(4, [-1.0 / sqrt(3.0),  1.0 / sqrt(3.0)], 1.0, deepcopy(material))
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
# 形状関数の定義
#----------------------------------------------------------------
function shape_functions(coordinate::Vector{Float64})

    xi::Float64, et::Float64 = coordinate

    N = [
        0.25 * (1.0 - xi) * (1.0 - et),
        0.25 * (1.0 + xi) * (1.0 - et),
        0.25 * (1.0 + xi) * (1.0 + et),
        0.25 * (1.0 - xi) * (1.0 + et)
    ]

    return N
end
#----------------------------------------------------------------
# ヤコビアン行列の計算
#----------------------------------------------------------------
function jacobian_matrix(nodes::Vector{Node}, coordinate::Vector{Float64})
    
    num_nodes::Int64 = length(nodes)

    # 形状関数の微分
    dN = shape_function_derivatives(coordinate)

    # 初期化
    J = zeros(2, 2)

    #for i in 1:num_nodes
    #    J += [nodes[element.nodes[i]].x, nodes[element.nodes[i]].y] * dN[i, :]
    #end

    return J
end
#----------------------------------------------------------------
# 形状関数の導関数の計算
#----------------------------------------------------------------
function shape_function_derivatives(coordinate::Vector{Float64})
    
    # 積分点座標の設定
    xi::Float64, et::Float64 = coordinate

    dN_dxi = [
        -0.25 * (1.0 - et) 
         0.25 * (1.0 - et) 
         0.25 * (1.0 + et) 
        -0.25 * (1.0 + et)
    ]

    dN_det = [
        -0.25 * (1.0 - xi) 
        -0.25 * (1.0 + xi) 
         0.25 * (1.0 + xi) 
         0.25 * (1.0 - xi)
    ]

    return [dN_dxi, dN_det]
end
#----------------------------------------------------------------
# 要素マトリクスの作成
#----------------------------------------------------------------
function compute_Ke(element)
    
    # 要素自由度
    num_dof = length(element.nodes)

    # 初期化
    Ke = zeros(Float64, num_dof, num_dof)

    # 積分点ループ
    for evaluate_point in element.evaluate_points

        # 形状関数の計算
        N = shape_functions(evaluate_point.coordinate)

        # ヤコビ行列の計算
        J = jacobian_matrix(element.nodes, evaluate_point.coordinate)

        # ガウス積分の重みとヤコビアンのスケーリング
        #weight = evaluate_point.weight * det(J)

        # ローカル剛性行列の計算
        #local_K += conductivity * (inv(J)' * inv(J)) * N * N' * weight

        # ローカル右辺ベクトルの計算（内部発熱項も含む）
        #local_F += N * element.heat_source * weight
    end

    return Ke
end
#----------------------------------------------------------------
# 要素マトリクスの作成
#----------------------------------------------------------------
function compute_Fbe(element)
    
    # 要素自由度
    num_equation = length(element.nodes)

    # 初期化
    Fe = zeros(Float64, num_equation, num_equation)

    # 積分点ループ
    for evaluate_point in element.evaluate_points

        # 形状関数の計算
        N = shape_functions(evaluate_point.coordinate)

        # ヤコビ行列の計算
        #J = jacobian_matrix(element.nodes, evaluate_point.coordinate)

        # ガウス積分の重みとヤコビアンのスケーリング
        #weight = gauss_weights[1] * gauss_weights[2] * det(J)

        # ローカル右辺ベクトルの計算（内部発熱項も含む）
        #local_F += N * element.heat_source * weight
    end

    return Fe
end