function hankel_transform(r²f::AbstractVector{T}, mesh::RadialMesh{T}, qs::AbstractVector{TT},
                          work_weights::AbstractVector{T}, work_integrand::AbstractVector{T},
                          l::Int=0;
                          quadrature_method::QuadratureMethod=Simpson(),
                          kwargs...)::AbstractVector{TT} where {T<:Real,TT<:Real}
    r = collect(mesh)
    jₗ = fast_sphericalbesselj(l)
    work_weights = @view work_weights[eachindex(r)]
    integration_weights!(work_weights, mesh, quadrature_method)
    work_integrand = @view work_integrand[eachindex(r)]

    f_q = map(qs) do q
        work_integrand .= r²f .* jₗ.(q .* r)
        integral = dot(work_weights, work_integrand)
        return 4TT(π) * integral
    end

    return f_q
end
