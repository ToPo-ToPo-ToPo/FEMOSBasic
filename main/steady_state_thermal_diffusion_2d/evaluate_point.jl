
include("material.jl")
#----------------------------------------------------------------
# 積分点構造体
#----------------------------------------------------------------
mutable struct EvaluatePoint
    id::Int64                       # 評価点の番号
    coordinate::Vector{Float64}     # 積分点座標
    weight::Float64                 # 積分点の重み
    material::Material              # 材料モデル
    #-------------------------------------------------------
    # 内部コンストラクタ
    #-------------------------------------------------------
    function EvaluatePoint(id::Int64, coordinate::Vector{Float64}, weight::Float64)
        return new(id, coordinate, weight, Material(0.0))
    end
    #-------------------------------------------------------
    # 内部コンストラクタ
    #-------------------------------------------------------
    EvaluatePoint(id::Int64, coordinate::Vector{Float64}, weight::Float64, material::Material) = new(id, coordinate, weight, material)
end