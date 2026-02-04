
"""
    EdgeMask{D}(A,f,delta)

Callable object that returns a edge mask function centered on
center of the domain, with `A`, the amplitude of the square wave,
`f` the frequency of the wave, and `delta` the with smoothing of the square wave

Examples
========

* Create a edge mask using a smoothed square wave, centered in the center of the domain.

```julia
julia> mask = EdgeMask{:z}(A, f, delta)
```
"""
struct EdgeMask{D, T}
     A :: T
     f :: T
     delta :: T
     Lx :: T
     Ly :: T
     threshold :: T

    function EdgeMask{D}(; A, f, delta, Lx, Ly, threshold) where D
        T = promote_type(typeof(A), typeof(f), typeof(delta), typeof(Lx), typeof(Ly), typeof(threshold) )
        return new{D, T}(A, f, delta, Lx, Ly, threshold)
    end
end    

@inline function gaussian(z,center, width)
    return exp(-(z - center)^2 / (2 * width^2))
end

@inline function (edge_mask::EdgeMask{:xy})(x, y, z)
    A = edge_mask.A
    δ = edge_mask.delta
    f = edge_mask.f
    Lx = edge_mask.Lx 
    Ly = edge_mask.Ly 
    threshold = edge_mask.threshold
 
    sin_x = sin(2π * A * f * x / (Lx/2) - π/2) / δ
    sin_y = sin(2π * A * f * y / (Ly/2) - π/2) / δ

    min_val = - 2 * atan( sin( -π/2 )/ δ )
    mask =  (atan(sin_x) + atan(sin_y)) + min_val

    normalized_mask = (2*mask/π)^2
    
    normalized_mask = normalized_mask #+ depth_mask

    if normalized_mask > 1 
	normalized_mask = 1
    end
    if normalized_mask < threshold
        normalized_mask = 0 
    end

    return normalized_mask
end
