
include("element.jl")
include("boundary_condition.jl")
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
        Ke = compute_Ke(element)

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
            #Fbe = compute_Fbe(element)

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