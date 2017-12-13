import COBRA, JuMP, GLPKMathProgInterface

export optimize


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