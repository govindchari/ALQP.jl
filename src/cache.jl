struct CACHE
    c_ineq::Array{Float64,1}
    c_ineq2::Array{Float64,1}
    c_eq::Array{Float64,1}
    c_n1::Array{Float64,1}
    c_n2::Array{Float64,1}
    c_n3::Array{Float64,1}
    c_ineq_n::Array{Float64,2}
    c_n_n::Array{Float64,2}

    function CACHE(n_eq, n_ineq, n)
        new(zeros(n_ineq), zeros(n_ineq), zeros(n_eq), zeros(n), zeros(n), zeros(n), zeros(n_ineq, n), zeros(n, n))
    end
end