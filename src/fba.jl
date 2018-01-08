import COBRA, JuMP, GLPKMathProgInterface

export optimize


"""
Optimizes a model. Returns status, f, v, where status
is the status of the LP optimization, f the objective 
value, and v the vector of fluxes.
"""
function optimize(model::COBRA.LPproblem)
    m,n = size(model.S)
    lp = JuMP.Model(solver=GLPKMathProgInterface.GLPKSolverLP())
    JuMP.@variable lp model.lb[k] ≤ x[k = 1:n] ≤ model.ub[k]

    for i = 1:m
        JuMP.@expression lp balance sum(model.S[i,k] * x[k] for k = 1:n)
        if model.csense[i] == '<' || model.csense[i] == 'L'
            JuMP.@constraint lp balance ≤ model.b[i]
        elseif model.csense[i] == '>' || model.csense[i] == 'G'
            JuMP.@constraint lp balance ≥ model.b[i]
        elseif model.csense[i] == '=' || model.csense[i] == 'E'
            JuMP.@constraint lp balance == model.b[i]
        else
            error("Invalid csense")
        end       
    end

    JuMP.@expression lp obj sum(model.c[k] * x[k] for k = 1:n)

    if model.osense == -1
        JuMP.@objective lp :Min obj
    elseif model.osense == 1
        JuMP.@objective lp :Max obj
    else
        error("Invalid osense")
    end

    sol = JuMP.solve(lp)

    return sol, JuMP.getobjectivevalue(lp), JuMP.getvalue(x)
end


"""
For each reaction in rxns, compute max. and min. flux attainable 
in the reaction. Returns two vectors fmin, fmax, where fmin 
contains the min fluxes and fmax the max fluxes.
"""
function fva(model::COBRA.LPproblem, rxns::AbstractVector{Int})
    m,n = size(model.S)
    @assert all(1 .≤ rxns .≤ n)

    lp = JuMP.Model(solver=GLPKMathProgInterface.GLPKSolverLP())
    JuMP.@variable lp model.lb[k] ≤ x[k = 1:n] ≤ model.ub[k]

    for i = 1:m
        JuMP.@expression lp balance sum(model.S[i,k] * x[k] for k = 1:n)
        if model.csense[i] == '<' || model.csense[i] == 'L'
            JuMP.@constraint lp balance ≤ model.b[i]
        elseif model.csense[i] == '>' || model.csense[i] == 'G'
            JuMP.@constraint lp balance ≥ model.b[i]
        elseif model.csense[i] == '=' || model.csense[i] == 'E'
            JuMP.@constraint lp balance == model.b[i]
        else
            error("Invalid csense")
        end       
    end

    fmin = zeros(length(rxns))
    fmax = zeros(length(rxns))

    for (i,k) in enumerate(rxns)
        JuMP.@objective lp :Min x[k]
        status = JuMP.solve(lp)
        @assert status == :Optimal
        fmin[i] = JuMP.getobjectivevalue(lp)

        JuMP.@objective lp :Max x[k]
        status = JuMP.solve(lp)
        @assert status == :Optimal
        fmax[i] = JuMP.getobjectivevalue(lp)
    end

    return fmin, fmax
end


"Perform fva for every reaction in the model."
fva(model::COBRA.LPproblem) = fva(model, 1 : length(model.rxns))

"Perform fva for the i'th reaction"
fva(model::COBRA.LPproblem, i::Int) = fva(model, [i])
