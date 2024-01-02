
using SparseArrays

#----------------------------------------------------------------
# ノード構造体
#----------------------------------------------------------------
struct Node
    x::Float64
    y::Float64
    temperature::Float64
    heat_source::Float64
end
#----------------------------------------------------------------
# 要素構造体
#----------------------------------------------------------------
struct Element
    nodes::Vector{Int}
    conductivity::Float64
end
#----------------------------------------------------------------
# 有限要素法で2次元熱伝導方程式を解く関数
#----------------------------------------------------------------
function solve_heat_conduction(nodes, elements, fixed_boundaries, heat_flux_boundaries)
    num_nodes = length(nodes)
    num_elements = length(elements)

    # グローバル剛性行列と右辺ベクトルの初期化
    K_global = sparse(zeros(Float64, num_nodes, num_nodes))
    F_global = zeros(Float64, num_nodes)

    # 要素ごとに剛性行列と右辺ベクトルを計算
    for elem_idx in 1:num_elements
        element = elements[elem_idx]
        local_K, local_F = compute_element_matrices(nodes, element)

        # グローバル剛性行列と右辺ベクトルに組み込み
        for i in 1:4
            for j in 1:4
                K_global[element.nodes[i], element.nodes[j]] += local_K[i, j]
            end
            F_global[element.nodes[i]] += local_F[i]
        end
    end

    # 温度固定境界条件の適用
    for (node_idx, temperature) in fixed_boundaries
        K_global[node_idx, :] .= 0.0
        K_global[node_idx, node_idx] = 1.0
        F_global[node_idx] = temperature
    end

    # 熱流束境界条件の適用
    for (node_idx, heat_flux) in heat_flux_boundaries
        F_global[node_idx] += heat_flux
    end

    # 疎行列を使って方程式を解く
    T = K_global \ F_global

    # 結果をノードに反映
    for i in 1:num_nodes
        nodes[i].temperature = T[i]
    end

    return nodes
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
    dN_dxi = shape_function_derivatives(xi, eta)

    J = zeros(2, 2)

    for i in 1:num_nodes
        J += [nodes[element.nodes[i]].x, nodes[element.nodes[i]].y] * dN_dxi[i, :]
    end

    return J
end
#----------------------------------------------------------------
# 形状関数の導関数の計算
#----------------------------------------------------------------
function shape_function_derivatives(xi, eta)
    dN_dxi = [-0.25 * (1 - eta), 0.25 * (1 - eta), 0.25 * (1 + eta), -0.25 * (1 + eta)]
    dN_deta = [-0.25 * (1 - xi), -0.25 * (1 + xi), 0.25 * (1 + xi), 0.25 * (1 - xi)]

    return [dN_dxi dN_deta]
end
#----------------------------------------------------------------
# 使用例
#----------------------------------------------------------------
# ノードの定義
nodes = [
    Node(0.0, 0.0, 0.0, 0.0),
    Node(1.0, 0.0, 0.0, 0.0),
    Node(1.0, 1.0, 0.0, 0.0),
    Node(0.0, 1.0, 0.0, 0.0)
]

# 要素の定義
elements = [
    Element([1, 2, 3, 4], 1.0)
]

# 境界条件の定義
fixed_boundaries = Dict(1 => 100.0, 2 => 0.0)
heat_flux_boundaries = Dict(3 => 10.0)

# 解の計算
result_nodes = solve_heat_conduction(nodes, elements, fixed_boundaries, heat_flux_boundaries)

# 結果の表示
for node in result_nodes
    println("Node $(node.x), $(node.y): Temperature = $(node.temperature)")
end
