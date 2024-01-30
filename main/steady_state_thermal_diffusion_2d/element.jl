
include("node.jl")
include("evaluate_point.jl")
using LinearAlgebra
#----------------------------------------------------------------
# 要素の構造体
#----------------------------------------------------------------
mutable struct Element
    id::Int64                                # 要素番号
    nodes::Vector{Node}                      # 要素節点
    thickness::Float64                       # 要素の厚み
    evaluate_points::Vector{EvaluatePoint}   # 積分点の情報

    #-------------------------------------------------------
    # 内部コンストラクタ
    #-------------------------------------------------------
    function Element(id::Int64, nodes::Vector{Node}, thickness::Float64, material::Material)
        
        # 評価点の作成
        evaluate_points = [
            EvaluatePoint(1, [-1.0 / sqrt(3.0), -1.0 / sqrt(3.0)], 1.0, deepcopy(material)),
            EvaluatePoint(2, [ 1.0 / sqrt(3.0), -1.0 / sqrt(3.0)], 1.0, deepcopy(material)),
            EvaluatePoint(3, [ 1.0 / sqrt(3.0),  1.0 / sqrt(3.0)], 1.0, deepcopy(material)),
            EvaluatePoint(4, [-1.0 / sqrt(3.0),  1.0 / sqrt(3.0)], 1.0, deepcopy(material))
        ]

        return new(id, nodes, thickness, evaluate_points)
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

    N1 = 0.25 * (1.0 - xi) * (1.0 - et)
    N2 = 0.25 * (1.0 + xi) * (1.0 - et)
    N3 = 0.25 * (1.0 + xi) * (1.0 + et)
    N4 = 0.25 * (1.0 - xi) * (1.0 + et)

    N = [N1 N2 N3 N4]

    return N
end
#----------------------------------------------------------------
# ヤコビアン行列の計算
#----------------------------------------------------------------
function make_jacobian(nodes::Vector{Node}, coordinate::Vector{Float64})
    
    # 要素節点数の取得
    num_node = length(nodes)

    # 形状関数の微分を計算
    dNdxi = shape_function_derivatives(coordinate)

    # 座標をまとめたマトリクスを作成
    coord_matrix = Matrix{Float64}(undef, num_node, 2)
    for i = 1 : num_node
        coord_matrix[i, :] = [nodes[i].coordinate[1] nodes[i].coordinate[2]]
    end

    # ヤコビ行列を計算する
    J = dNdxi * coord_matrix

    return J
end
#----------------------------------------------------------------
# 形状関数の導関数の計算
#----------------------------------------------------------------
function shape_function_derivatives(coordinate::Vector{Float64})
    
    # 積分点座標の設定
    xi, et = coordinate

    # 形状関数の微分を計算
    dNdxi = Matrix{Float64}(undef, 2, 4)

    dNdxi[1, :] = [-0.25 * (1.0 - et) 0.25 * (1.0 - et) 0.25 * (1.0 + et) -0.25 * (1.0 + et)]
    dNdxi[2, :] = [-0.25 * (1.0 - xi) -0.25 * (1.0 + xi) 0.25 * (1.0 + xi) 0.25 * (1.0 - xi)]

    return dNdxi
end
#----------------------------------------------------------------
# Bマトリクスの作成
#----------------------------------------------------------------
function make_B(nodes, coordinate)

    # 形状関数の微分を計算
    dNdxi = shape_function_derivatives(coordinate)

    # ヤコビ行列の計算
    J = make_jacobian(nodes, coordinate)

    # dNdxの計算: dNdxi = J * dNdx -> dNdx = J_invers * dNdxi
    dNdx = J \ dNdxi

    # Bマトリクスの整形
    B = dNdx

    return B

end
#----------------------------------------------------------------
# 要素マトリクスの作成
#----------------------------------------------------------------
function make_Ke(element)
    
    # 要素自由度
    num_dof = length(element.nodes)

    # 初期化
    Ke = zeros(Float64, num_dof, num_dof)

    # 積分点ループ
    for evaluate_point in element.evaluate_points

        # ヤコビ行列の計算
        J = make_jacobian(element.nodes, evaluate_point.coordinate)

        # Bマトリクスの計算
        B = make_B(element.nodes, evaluate_point.coordinate)

        # ローカル剛性行列の計算
        Ke += evaluate_point.material.conductivity * B' * B * evaluate_point.weight * det(J) * element.thickness

    end

    return Ke
end
#----------------------------------------------------------------
# 要素のソース項の作成
#----------------------------------------------------------------
function make_Fbe(element, value)
    
    # 要素自由度
    num_dof = length(element.nodes)

    # 初期化
    Fe = zeros(Float64, num_dof)

    # 積分点ループ
    for evaluate_point in element.evaluate_points

        # 形状関数の計算
        N = shape_functions(evaluate_point.coordinate)

        # ヤコビ行列の計算
        J = make_jacobian(element.nodes, evaluate_point.coordinate)

        # 要素ベクトルの計算
        Fe += N' * value * evaluate_point.weight * det(J) * element.thickness
    
    end

    return Fe
end