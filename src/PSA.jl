module PSA

using DataFrames, CSV, LightGraphs

export Network, import_network, idx, rev_idx, select_names, select_by, idx_by, to_symbol, append_idx_col!

# include("auxilliaries.jl") already in lopf
include("lopf.jl")

mutable struct network_mutable
    buses::DataFrame
    generators::DataFrame
    loads::DataFrame
    lines::DataFrame
    links::DataFrame
    storage_units::DataFrame
    stores::DataFrame
    transformers::DataFrame
    carriers::DataFrame
    global_constraints::DataFrame
    buses_t::Dict{String,DataFrame}
    generators_t::Dict{String,DataFrame}
    loads_t::Dict{String,DataFrame}
    lines_t::Dict{String,DataFrame}
    links_t::Dict{String,DataFrame}
    storage_units_t::Dict{String,DataFrame}
    stores_t::Dict{String,DataFrame}
    transformers_t::Dict{String,DataFrame}
    snapshots::DataFrame
end



function Network(
# static
    buses = DataFrame(repeat([Bool[]], outer=16),
            [:name, :v_nom, :type, :x, :y, :carrier, :v_mag_pu_set, :v_mag_pu_min,
            :v_mag_pu_max, :control, :sub_network, :p, :q,
            :v_mag_pu, :v_ang, :marginal_price]),

    carriers = DataFrame(repeat([Bool[]], outer=2),
                 [:name, :co2_emissions]),

    generators = DataFrame(repeat([Bool[]], outer=31),
            [:name, :bus, :control, :type, :p_nom, :p_nom_extendable,
            :p_nom_min, :p_nom_max, :p_min_pu, :p_max_pu, :p_set, :q_set,
            :sign, :carrier, :marginal_cost, :capital_cost, :efficiency,
            :committable, :start_up_cost, :shut_down_cost, :min_up_time,
            :min_down_time, :initial_status, :ramp_limit_up, :ramp_limit_down,
            :ramp_limit_start_up, :ramp_limit_shut_down, :p, :q, :p_nom_opt, :status]),

    global_constraints = DataFrame(repeat([Bool[]], outer=6),
            [:name, :type, :carrier_attribute, :sense, :constant, :mu]),

    line_types = DataFrame(repeat([Bool[]], outer=9),
            [:name, :f_nom, :r_per_length, :x_per_length, :c_per_length,
            :i_nom, :mounting, :cross_section, :references]),

    lines = DataFrame(repeat([Bool[]], outer=33),
            [:name, :bus0, :bus1, :type, :x, :r, :g, :b,
            :s_nom, :s_nom_extendable, :s_nom_min, :s_nom_max,
            :s_max_pu, :capital_cost, :length, :terrain_factor, :num_parallel,
            :v_ang_min, :v_ang_max, :sub_network, :p0, :q0, :p1, :q1,
            :x_pu, :r_pu, :g_pu, :b_pu, :x_pu_eff, :r_pu_eff,
            :s_nom_opt, :mu_lower, :mu_upper]),

    links = DataFrame(repeat([Bool[]], outer=21),
            [:name, :bus0, :bus1, :type, :efficiency, :p_nom, :p_nom_extendable,
            :p_nom_min, :p_nom_max, :p_set, :p_min_pu, :p_max_pu, :capital_cost,
            :marginal_cost, :length, :terrain_factor, :p0, :p1,
            :p_nom_opt, :mu_lower, :mu_upper]),

    loads = DataFrame(repeat([Bool[]], outer=8),
            [:name, :bus, :type, :p_set, :q_set, :sign, :p, :q]),

    shunt_impendances = DataFrame(repeat([Bool[]], outer=9),
            [:name, :bus, :g, :b, :sign, :p, :q, :g_pu, :b_pu]),

    storage_units = DataFrame(repeat([Bool[]], outer=29),
            [:name, :bus, :control, :type, :p_nom, :p_nom_extendable, :p_nom_min,
            :p_nom_max, :p_min_pu, :p_max_pu, :p_set, :q_set, :sign, :carrier,
            :marginal_cost, :capital_cost, :state_of_charge_initial, :state_of_charge_set,
            :cyclic_state_of_charge, :max_hours, :efficiency_store, :efficiency_dispatch,
            :standing_loss, :inflow, :p, :q, :state_of_charge, :spill, :p_nom_opt]),

    stores = DataFrame(repeat([Bool[]], outer=22),
            [:name, :bus, :type, :e_nom, :e_nom_extendable, :e_nom_min,
            :e_nom_max, :e_min_pu, :e_max_pu, :e_initial, :e_cyclic, :p_set,
            :cyclic_state_of_charge, :q_set, :sign, :marginal_cost,
            :capital_cost, :standing_loss, :p, :q, :e, :e_nom_opt]),

    transformer_types = DataFrame(repeat([Bool[]], outer=16),
            [:name, :f_nom, :s_nom, :v_nom_0, :v_nom_1, :vsc, :vscr, :pfe,
            :i0, :phase_shift, :tap_side, :tap_neutral, :tap_min, :tap_max, :tap_step, :references]),

    transformers = DataFrame(repeat([Bool[]], outer=36),
            [:name, :bus0, :bus1, :type, :model, :x, :r, :g, :b, :s_nom,
            :s_nom_extendable, :s_nom_min, :s_nom_max, :s_max_pu, :capital_cost,
            :num_parallel, :tap_ratio, :tap_side, :tap_position, :phase_shift,
            :v_ang_min, :v_ang_max, :sub_network, :p0, :q0, :p1, :q1, :x_pu,
            :r_pu, :g_pu, :b_pu, :x_pu_eff, :r_pu_eff, :s_nom_opt, :mu_lower, :mu_upper]),

# time_dependent
    buses_t=Dict{String,DataFrame}([("marginal_price",DataFrame()), ("v_ang", DataFrame()),
            ("v_mag_pu_set",DataFrame()), ("q", DataFrame()),
            ("v_mag_pu", DataFrame()), ("p", DataFrame()), ("p_max_pu", DataFrame())]),
    generators_t=Dict{String,DataFrame}(
            [("marginal_price",DataFrame()), ("v_ang", DataFrame()),
            ("v_mag_pu_set",DataFrame()), ("q", DataFrame()),
            ("v_mag_pu", DataFrame()), ("p", DataFrame()),
            ("p_min_pu", DataFrame()), ("p_max_pu", DataFrame())
            ]),
    loads_t=Dict{String,DataFrame}([
            ("q_set",DataFrame()), ("p_set", DataFrame()),
            ("q", DataFrame()), ("p", DataFrame())
            ]),
    lines_t=Dict{String,DataFrame}([("q0",DataFrame()), ("q1", DataFrame()),
            ("p0",DataFrame()), ("p1", DataFrame()),
            ("mu_lower", DataFrame()), ("mu_upper", DataFrame())]),
    links_t=Dict{String,DataFrame}([("p_min_pu",DataFrame()), ("p_max_pu", DataFrame()),
            ("p0",DataFrame()), ("p1", DataFrame()),
            ("mu_lower", DataFrame()), ("mu_upper", DataFrame()),
            ("efficiency",DataFrame()), ("p_set", DataFrame())]),
    storage_units_t=Dict{String,DataFrame}([("p_min_pu",DataFrame()), ("p_max_pu", DataFrame()),
        ("inflow",DataFrame()), ("mu_lower", DataFrame()), ("mu_upper", DataFrame()),
        ("efficiency",DataFrame()), ("p_set", DataFrame())]),
    stores_t=Dict{String,DataFrame}([("e_min_pu",DataFrame()), ("e_max_pu", DataFrame()),
        ("inflow",DataFrame()), ("mu_lower", DataFrame()), ("mu_upper", DataFrame()),
        ("efficiency",DataFrame()), ("e_set", DataFrame())]),
    transformers_t= Dict{String,DataFrame}( [("q1", DataFrame()), ("q0", DataFrame()), ("p0", DataFrame()),
        ("p1", DataFrame()), ("mu_upper", DataFrame()), ("mu_lower", DataFrame()),
        ("s_max_pu", DataFrame())]),
    snapshots=DataFrame([Bool[]], [:t])

    )
    network_mutable(
        buses, generators, loads, lines, links, storage_units, stores, transformers, carriers,
        global_constraints,
        buses_t, generators_t, loads_t, lines_t, links_t, storage_units_t, stores_t, transformers_t,
        snapshots)
end


function import_network(folder)
    network = Network()
    !ispath("$folder") ? error("Path not existent") : nothing
    components = static_components(network)
    for component=components
        componentpath = joinpath("$folder", "$component.csv")
        if ispath(componentpath)
            # fallback for missing values
            try
                setfield!(network,component,CSV.read(componentpath, DataFrame; truestrings=["True"],
                                                    falsestrings=["False"]))
            catch y
                if (typeof(y)==Missings.MissingException) | (typeof(y) == BoundsError)
                    setfield!(network,component,readtable(componentpath; truestrings=["True"],
                                                    falsestrings=["False"]))

                end
            end
        end
        # convert any ints to floats
        df = getfield(network, component)
        for name in names(df)
            if name == :name
                nothing
            elseif typeof(df[!, name]) == Array{Union{Int64, Missings.Missing},1}
                println("converting column $name of $component from Int to Float")
                df[name] = float(df[name])
            end
        end
    end
    components_t = time_dependent_components(network)
    for component_t=components_t
        for attr in keys(getfield(network, component_t))
            component = Symbol(String(component_t)[1:end-2])
            if ispath("$folder/$component-$attr.csv")
                # fallback for missing values for a non-null column type, might be deprecated soon
                try
                    getfield(network,component_t)[attr]= (
                    CSV.read("$folder/$component-$attr.csv"; truestrings=["True"], falsestrings=["False"]) )
                catch y
                    if (typeof(y)==Missings.MissingException) | (typeof(y) == BoundsError)
                        getfield(network,component_t)[attr]= (
                            readtable("$folder/$component-$attr.csv"; truestrings=["True"], falsestrings=["False"]) )
                    end
                end
            end
        end
    end
    initializer = Network()
    for field=setdiff(fieldnames(typeof(network)), fieldnames(typeof(initializer)))
        setfield!(network, field, getfield(initializer, field))
    end
    return network
end


function export_network(network, folder)
    !ispath(folder) ? mkdir(folder) : nothing
    components_t = [field for field=fieldnames(network) if String(field)[end-1:end]=="_t"]
    for field_t in components_t
        for df_name=keys(getfield(network,field_t))
            if nrow(getfield(network,field_t)[df_name])>0
                field = Symbol(String(field_t)[1:end-2])
                writetable("$folder/$field-$df_name.csv", getfield(network,field_t)[df_name])
            end
        end
    end
    components = [field for field=fieldnames(network) if String(field)[end-1:end]!="_t"]
    for field in components
        writetable("$folder/$field.csv", getfield(network, field))
    end
end



end
