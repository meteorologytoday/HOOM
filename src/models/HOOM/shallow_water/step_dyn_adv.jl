
function advectDynamic!(
    model   :: Model,
    Δt      :: Float64,
)

    # 1. derive barotropic and baroclinic flow
    # 2. derive G of each terms
    #
    #    a)
    #
    state = model.state
    tcr_adv = model.tcr_adv
    env = model.env

    cal∇b!
    calCof!
    calv∇v!
    calG_bc!

end

