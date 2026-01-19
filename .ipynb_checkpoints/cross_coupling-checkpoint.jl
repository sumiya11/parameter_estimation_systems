let

# Independent variable:
@parameters t

# Parameters:
ps = @parameters k₁ k₂ k₃ Initial(S₁C(t)) Initial(Cˍt(t)) Initial(S₂ˍt(t)) Initial(S₁Cˍt(t)) Initial(S₂(t)) Initial(P(t)) Initial(Pˍt(t)) Initial(S₁ˍt(t)) Initial(C(t)) Initial(CP(t)) Initial(S₁(t)) Initial(CPˍt(t))

# Species:
sps = @species S₁(t) C(t) S₁C(t) S₂(t) CP(t) P(t)

# Reactions:
rxs = [
	Reaction(k₁, [S₁, C], [S₁C], [1, 1], [1]),
	Reaction(k₂, [S₁C, S₂], [CP], [1, 1], [1]),
	Reaction(k₃, [CP], [C, P], [1], [1, 1])
]

# Defaults:
defaults = Dict([Initial(S₂) => false, Initial(C) => false, Initial(S₁C) => false, Initial(Cˍt) => false, Initial(CPˍt) => false, Initial(S₁ˍt) => false, Initial(P) => false, Initial(S₁) => false, Initial(S₁Cˍt) => false, Initial(Pˍt) => false, Initial(S₂ˍt) => false, Initial(CP) => false])

# Declares ReactionSystem model:
rs = ReactionSystem(rxs, t, sps, ps; name = Symbol("##ReactionSystem#371"), defaults)
complete(rs)

end