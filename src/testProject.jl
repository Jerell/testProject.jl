module testProject
using Unitful

struct Pipe
  length::typeof(1u"m")
  diameter::typeof(1u"m")
  height::typeof(1u"m")
end

struct Fluid
  pressure::typeof(1u"bar")
  temperature::typeof(1u"°C")
  flowrate::typeof(1u"kg/s")
end

struct EoSParams{T<:Real}
  compound::String
  b::Quantity{T,dimension(u"L/mol"),typeof(u"L/mol")}
  a₀::Quantity{T,dimension(u"bar * L^2 *mol^-2"),typeof(u"bar * L^2 *mol^-2")}
  c₁::T
  ϵ::Quantity{T,dimension(u"bar*L*mol^-1"),typeof(u"bar*L*mol^-1")}
  β::T
  Tᶜᵐ::Quantity{T,dimension(u"K"),typeof(u"K")}
  Pᶜᵐ::Quantity{T,dimension(u"bar"),typeof(u"bar")}
  mₘ::T

  EoSParams(name::String, b::Real, a₀::Real, c₁::Real, ϵ::Real, β::Real, Tᶜᵐ::Real, Pᶜᵐ::Real, mₘ::Real) = new{typeof(b)}(
    name,
    b * u"L/mol",
    a₀ * u"bar * L^2 *mol^-2",
    c₁,
    ϵ * u"bar*L*mol^-1",
    β * 10^3,
    Tᶜᵐ * u"K",
    Pᶜᵐ * u"bar",
    mₘ,
  )
end


function α(params::EoSParams)
  return function (T::Quantity{Real,dimension(u"K"),typeof(u"K")})
    params.a₀ * (1 + params.c₁ * (1 - √(T / params.Tᶜᵐ)))^2
  end
end


Ωₐ = 0.42748
Ωᵦ = 0.08664

compoundparams = [
  EoSParams("methanol", 0.0309, 4.0531, 0.4310, 245.91, 16.1, 363.75, 84.587, 0.33997)
  EoSParams("ethanol", 0.0491, 8.6716, 0.7369, 215.32, 8.0, 463.34, 68.005, 0.67459)
  EoSParams("propanol", 0.0641, 11.9102, 0.9171, 210.00, 8.1, 490.26, 55.070, 0.84223)
  EoSParams("2-propanol", 0.0641, 10.6019, 0.9468, 210.00, 9.1, 449.76, 50.544, 0.84323)
  EoSParams("butanol", 0.0797, 15.6949, 0.9784, 210.00, 8.2, 518.61, 46.874, 0.90340)
  EoSParams("pentanol", 0.0974, 22.7576, 0.9358, 210.00, 3.6, 577.32, 42.674, 0.92220)
  EoSParams("octanol", 0.1485, 41.5822, 1.1486, 267.59, 0.14, 666.35, 32.335, 1.17497)
  EoSParams("formic acid", 0.0300, 5.6336, 0.3338, 419.17, 15.5, 485.96, 116.688, 0.29452)
  EoSParams("acetic acid", 0.0468, 9.1196, 0.4644, 403.23, 4.5, 508.03, 78.198, 0.41601)
  EoSParams("propanoic acid", 0.0641, 13.2676, 0.6891, 399.75, 2.1, 540.81, 60.777, 0.63149)
  EoSParams("(mono)ethylene glycol, MEG", 0.0514, 10.8190, 0.6744, 197.52, 14.1, 584.11, 81.862, 0.56931)
  EoSParams("diethylene glycol, DEG", 0.0921, 26.4080, 0.7991, 196.84, 6.4, 718.69, 56.213, 0.77422)
  EoSParams("triethylene glycol, TEG", 0.1321, 39.1260, 1.1692, 143.37, 18.8, 747.04, 40.737, 1.13254)
  EoSParams("propylene glycol, PG", 0.0675, 13.8360, 0.9372, 174.42, 19.0, 555.47, 59.280, 0.83731)
  EoSParams("methylamine", 0.0361, 5.4886, 0.6044, 114.66, 33.7, 391.47, 78.116, 0.56108)
  EoSParams("ethylamine", 0.0531, 9.0567, 0.7382, 93.175, 43.0, 432.23, 58.637, 0.70476)
  EoSParams("diethylamine", 0.0857, 17.2460, 0.8838, 37.008, 110.8, 493.37, 41.471, 0.87840)
  EoSParams("water", 0.0145, 1.2277, 0.6736, 166.55, 69.2, 303.17, 150.459, 0.38016)
]

function ŋ(params::EoSParams)
  return function (Vₘ::Quantity{Real,dimension(u"L/mol"),typeof(u"L/mol")})
    return params.b / (4Vₘ)
  end
end

function g(params::EoSParams)
  return function (Vₘ::Quantity{Real,dimension(u"L/mol"),typeof(u"L/mol")})
    return 1 / (1 - 1.9ŋ(params)(Vₘ))
  end
end


function eos(params::EoSParams)
  dlngdvm = 1
  R = 0.0831446261815324
  return (
    T::Quantity{Real,dimension(u"K"),typeof(u"K")},
    Vₘ::Quantity{Real,dimension(u"L/mol"),typeof(u"L/mol")}
  ) -> begin
    Δᴬᴮ = g(params)(Vₘ) *
          (exp(params.ϵ / (R * T)) - 1) *
          params.b * params.β

    return R * T / (Vₘ - params.b) -
           α(params)(T) / (Vₘ * (Vₘ + params.b)) -
           1 / 2 * (R * T / Vₘ^-1) * (1 + Vₘ^-1 * dlngdvm)
    # * association term
  end
end

export Pipe, Fluid, EoSParams, compoundparams

end
