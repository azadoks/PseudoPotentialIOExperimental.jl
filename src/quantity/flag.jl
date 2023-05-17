abstract type PsPQuantityFlag end
Base.Broadcast.broadcastable(flag::PsPQuantityFlag) = Ref(flag)

struct LocalPotential <: PsPQuantityFlag end

abstract type ProjectorFlag <: PsPQuantityFlag end
struct BetaProjector <: ProjectorFlag end
struct ChiProjector <: ProjectorFlag end

struct BetaCoupling <: PsPQuantityFlag end

struct ValenceChargeDensity <: PsPQuantityFlag end
struct CoreChargeDensity <: PsPQuantityFlag end

struct AugmentationFunction <: PsPQuantityFlag end
struct AugmentationCoupling <: PsPQuantityFlag end
