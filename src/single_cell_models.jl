function generator_OU(mother_x,theta)
    # simple OU model with jumps at division 
    # time, expr., growth_rate, size
    τ,D,α,σY,dt = theta
    # time, growth rate,size
    X = [[mother_x[1]+dt,mother_x[2],mother_x[3]/2]] 
    target_fold_change = log(2) - α*(log(X[end][3])) + rand(Normal(0,σY))
    while X[end][3]< exp(target_fold_change)*X[1][3] 
        X_new = zeros(3)
        time,gr,l = copy(X[end][:])
        X_new[1] = time+dt
        X_new[2] = gr + 1/τ*(1-gr)*dt +  sqrt(2*D*dt)*randn() 
        X_new[3] = l + gr*l*dt
        push!(X,X_new)
    end
    return Matrix(hcat(X...)')
end

function generator_OU_rate(mother_x,theta)
    # simple OU model with jumps at division 
    # time, expr., growth_rate, size
    τ,D,β,dt = theta
    # time, growth rate,size
    X = [[mother_x[1]+dt,mother_x[2],mother_x[3]/2]] 
    y0 = log.(mother_x[3]/2)
    #target_fold_change = log(2) - α*(log(X[end][3])) + rand(Normal(0,σY))
    div = false
    while !div 
        X_new = zeros(3)
        time,gr,l = copy(X[end][:])
        X_new[1] = time+dt
        X_new[2] = gr + 1/τ*(1-gr)*dt +  sqrt(2*D*dt)*randn() 
        X_new[3] = l + max(10e-5,gr)*l*dt
        push!(X,X_new)
        div_prob = β(gr,log.(l),y0)*dt
        r = rand(Uniform(0,1))
        if r<div_prob
            div = true
        end
    end
    return Matrix(hcat(X...)')
end


function generator_RG(mother_x,theta)
    # random growth rate model
    # time, expr., growth_rate, size
    sgr,α,σY,dt,a = theta
    # time, growth rate,size
    gr = rand(Normal(1 .- mother_x[3]/2*a,sgr))
    gr = max(0.1,gr)
    X = [[mother_x[1]+dt,gr,mother_x[3]/2]] 
    target_fold_change = log(2) - α*(log(X[end][3])) + rand(Normal(0,σY))
    while X[end][3]< exp(target_fold_change)*X[1][3] 
        X_new = zeros(3)
        time,gr,l = copy(X[end][:])
        X_new[1] = time+dt
        X_new[2] = gr 
        X_new[3] = l + gr*l*dt
        push!(X,X_new)
    end
    return Matrix(hcat(X...)')
end

function φ_gaussian(y, y0, σY, α)
    μ = log(2) - α * y0
    A = exp(-((y - μ)^2) / (2 * σY^2))/(2 * π * σY^2)^(1/2)
    B = (1 - erf((y - μ) / (sqrt(2) * σY)))
    return A / B
end
