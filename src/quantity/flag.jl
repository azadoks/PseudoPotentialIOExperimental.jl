abstract type AtomicQuantityFlag end
Base.Broadcast.broadcastable(flag::AtomicQuantityFlag) = Ref(flag)

struct LocalPotential <: AtomicQuantityFlag end

abstract type ProjectorFlag <: AtomicQuantityFlag end
struct NonLocalProjector <: ProjectorFlag end
abstract type StateFlag <: ProjectorFlag end
struct PseudoState <: StateFlag end
struct HydrogenicState <: StateFlag end

struct NonLocalCoupling <: AtomicQuantityFlag end

abstract type DensityFlag <: AtomicQuantityFlag end
struct CoreDensity <: DensityFlag end
abstract type ValenceDensityFlag <: DensityFlag end
struct PseudoValenceDensity <: ValenceDensityFlag end
struct GaussianValenceDensity <: ValenceDensityFlag end

struct AugmentationFunction <: AtomicQuantityFlag end
struct AugmentationCoupling <: AtomicQuantityFlag end

# Quantity
#  - Local potential
#  - Projector
#    - Non-local potential projector
#    - State projector
#      - Pseudo-atomic state
#      - Hydrogenic state
#  - Charge density
#    - Core charge density
#    - Valence charge density
#      - Pseudo-atomic valence charge density
#      - Gaussian valence charge density
#  - Augmentation function
