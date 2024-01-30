
include("node.jl")
include("evaluate_point.jl")
#----------------------------------------------------------------
# 線要素の構造体
#----------------------------------------------------------------
mutable struct Line
    id::Int64                                # 要素番号
    nodes::Vector{Node}                      # 要素節点
    cross_section::Float64                   # 要素の太さ
end
#----------------------------------------------------------------
# 要素長さを計算
#----------------------------------------------------------------
function compute_length(line::Line)
    norm = 0.0
    for i = 1 : 3
        norm += (line.nodes[2].coordinate[i] - line.nodes[1].coordinate[i])^2.0
    end
    return sqrt(norm)
end
#----------------------------------------------------------------
# 要素中心座標を計算する関数
#----------------------------------------------------------------
function compute_element_center(element::Line)

    # 初期化
    center_coord = zeros(Float64, 3)

    # 主な処理
    for node in element.nodes
        center_coord += node.coordinate
    end

    # 平均化した値を返す
    return center_coord / length(element.nodes)
end