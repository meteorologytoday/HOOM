

module HOOM

    dynamic_core    # Non parallizable
    tracer_core     # parallizable
    mld_core        # parallizable

    data_manager    # manage underlying data exchange
    data_recorder   # history file tool

end
