
include("node.jl")
include("evaluate_point.jl")
#----------------------------------------------------------------
# 線要素の構造体
#----------------------------------------------------------------
mutable struct Line
    id::Int64                                # 要素番号
    nodes::Vector{Node}                      # 要素節点
    cross_section::Float64                   # 要素の太さ
    evaluate_points::Vector{EvaluatePoint}   # 積分点の情報
    #-------------------------------------------------------
    # 内部コンストラクタ
    #-------------------------------------------------------
    function Line(id::Int64, nodes::Vector{Node}, cross_section::Float64)
        
        # 評価点の作成
        evaluate_points = [
            EvaluatePoint(1, [-1.0 / sqrt(3.0)], 1.0),
            EvaluatePoint(2, [ 1.0 / sqrt(3.0)], 1.0)
        ]

        return new(id, nodes, cross_section, evaluate_points)
    end
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