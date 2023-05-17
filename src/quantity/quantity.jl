abstract type AbstractPsPQuantity end
Base.Broadcast.broadcastable(quantity::AbstractPsPQuantity) = Ref(quantity)
