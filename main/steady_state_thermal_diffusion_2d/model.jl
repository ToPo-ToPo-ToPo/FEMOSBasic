
include("element.jl")
include("line.jl")
include("boundary_condition.jl")
using SparseArrays
#----------------------------------------------------------------
# モデル
#----------------------------------------------------------------
mutable struct Model
    nodes::Vector{Node}                       # 全節点の情報
    elements::Vector{Element}                 # 全要素の情報
    lines::Vector{Line}                       # 線要素の情報
end
#----------------------------------------------------------------
# 有限要素法で2次元熱伝導方程式を解く関数
#----------------------------------------------------------------
function solve(model::Model, dirichlet_bcs::Vector{Dirichlet}, neumman_bcs::Vector{Neumman}, sources::Vector{Source})
    
    # 基本データの取得
    nodes = model.nodes
    elements = model.elements
    
    # 総自由度の取得
    num_total_dof = length(nodes)

    # グローバル行列と全体ベクトルの初期化
    K  = spzeros(Float64, num_total_dof, num_total_dof)
    Fq = zeros(Float64, num_total_dof)
    Fb = zeros(Float64, num_total_dof)

    # 全体熱流束ベクトルを作成
    for neumman_bc in neumman_bcs
        # 熱流束の値を取得
        value = neumman_bc.value
        
        # ベクトルに足し込む
        for line in neumman_bc.lines
            # 節点に与える熱流束
            heat_flux = 0.5 * value * line.cross_section * compute_length(line)
            
            # アセンブリング
            for i = 1 : length(line.nodes)
                # 全体系での自由度番号を取得
                idof = line.nodes[i].id
                Fq[idof] += heat_flux
            end
        end
    end

    # 全体熱伝導マトリクスを作成
    for element in elements
        # 要素行列作成
        Ke = make_Ke(element)

        # アセンブリング
        for i = 1 : length(element.nodes)
            
            # 全体系での自由度番号を取得
            idof = element.nodes[i].id
            
            for j = 1 : length(element.nodes)
                
                # 全体系での自由度番号を取得
                jdof = element.nodes[j].id
                
                # 全体マトリクスに足し込み
                K[idof, jdof] += Ke[i, j]
            end
        end
    end

    # ソース項ベクトルを作成
    for source in sources
        value = source.value
        for element in source.elements
            # 要素ベクトルを計算
            Fbe = make_Fbe(element, value)

            # アセンブリング
            for i = 1 : length(element.nodes)
                # 全体系での自由度番号を取得
                idof = element.nodes[i].id

                # 足し込み
                Fb[idof] += Fbe[i]
            end
        end
    end

    # 方程式の作成
    lhs = K
    rhs = Fq + Fb

    # ディレクレ境界条件の適用
    for dirichlet_bc in dirichlet_bcs
        # 条件値を取得
        value = dirichlet_bc.value
        
        # 条件が与えられる節点ループ
        for node in dirichlet_bc.nodes
            # 条件を与える自由度番号を取得
            idof = node.id

            # value != 0の場合の処置
            for jdof = 1 : num_total_dof
                rhs[jdof] -= value * lhs[jdof, idof]
            end

            # 左辺マトリクスの修正
            lhs[:,    idof] = zeros(num_total_dof)
            lhs[idof,    :] = zeros(num_total_dof)
            lhs[idof, idof] = 1.0

            # 右辺ベクトルの修正
            rhs[idof] = value
        end
    end
    
    # デフォルトのスパースソルバーの場合
    T = lhs \ rhs

    return T
end