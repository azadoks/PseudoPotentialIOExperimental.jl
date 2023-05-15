abstract type PsPQuantityFlag end

struct LocalPotential <: PsPQuantityFlag end

abstract type ProjectorFlag <: PsPQuantityFlag end
struct BetaProjector <: ProjectorFlag end
struct ChiProjector <: ProjectorFlag end

struct ValenceChargeDensity <: PsPQuantityFlag end
struct CoreChargeDensity <: PsPQuantityFlag end

struct AugmentationFunction <: PsPQuantityFlag end
