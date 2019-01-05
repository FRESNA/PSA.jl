using JuMP
using MathProgBase

include("utils.jl")

function build_lopf(network, solver; rescaling::Bool=false,formulation::String="angles_linear",
                    investment_type::String="continuous",
                    blockmodel::Bool=false, benders::String="", snapshot_number=0, N_groups=1, blockstructure::Bool=false)

    # This function is organized as the following:
    #
    # 0.        Initialize model
    # 1. - 5.   add generators,  lines, links, storage_units,
    #           stores to the model:
    #               .1 separate different types from each other
    #               .2 define number of different types
    #               .3 add variables to the model
    #               .4 set contraints for extendables
    #               .5 set charging constraints (storage_units and stores)
    # 6.        give power flow formulation
    # 7.        set global constraints
    # 8.        give objective function
    
    # Conventions:
    #   - all variable names start with capital letters
    #   - Additional capital words - which are not variables - are N (number of buses),
    #     T (number of snapshots), L (number of lines) defining the number of considered variables added to the model.
    
    # sanity checks
    snapshot_number>0 && benders!="slave" ? error("Can only specify one single snapshot for slave-subproblem!") : nothing
    blockmodel && benders!="" ? error("Can either do manual benders decomposition or use BlockDecomposition.jl!") : nothing
    
    # TODO: temporarily disable rescaling
    if rescaling
        rescaling = false
        println("Disabled rescaling contrary to setting!")
    end

    rf = 1
    rescaling_dict = rescaling_factors(rescaling)

    N = nrow(network.buses)
    L = nrow(network.lines)
    T = nrow(network.snapshots) #normally snapshots
    calculate_dependent_values!(network)
    nrow(network.loads_t["p"])!=T ? network.loads_t["p"]=network.loads_t["p_set"] : nothing
    candidates = line_extensions_candidates(network)

    # shortcuts
    buses = network.buses
    generators = network.generators
    links = network.links
    storage_units = network.storage_units
    stores = network.stores
    lines = network.lines

    # indices
    reverse_busidx = rev_idx(buses)
    busidx = idx(buses)
    reverse_lineidx = rev_idx(lines)
    lineidx = idx(lines)

    # # iterator bounds
    N_fix_G = sum(((.!generators[:p_nom_extendable]) .& (.!generators[:commitable])))
    N_ext_G = sum(convert(BitArray, generators[:p_nom_extendable]))
    N_com_G = sum(convert(BitArray, generators[:commitable]))
    N_fix_LN = sum((.!lines[:s_nom_extendable]))
    N_ext_LN = sum(.!(.!lines[:s_nom_extendable]))
    N_fix_LK = sum(.!links[:p_nom_extendable])
    N_ext_LK = sum(.!(.!links[:p_nom_extendable]))
    N_fix_SU = sum(.!storage_units[:p_nom_extendable])
    N_ext_SU = sum(.!(.!storage_units[:p_nom_extendable]))
    N_SU = nrow(storage_units)
    N_fix_ST = sum(.!stores[:e_nom_extendable])
    N_ext_ST = sum(.!(.!stores[:e_nom_extendable]))
    N_ST = N_fix_ST + N_ext_ST
    
    if N_ext_LN > 0
        error("Code currently only works if all lines are extendable lines!")
    end
    
    if N_com_G > 0
        println("WARNING, no unit commitment yet")
    end

    sn = snapshot_number
    if sn>0
        Te = 1 
        Trange = sn:sn
    else
        Te = T
        Trange = 1:T
    end

    nt=length(Trange)
    
    if blockmodel
        m = BlockModel(solver=solver)
    else
        m = Model(solver=solver)
    end

# --------------------------------------------------------------------------------------------------------

# 0. add all variables to the Model

    if benders != "slave"
        @variable(m, G_p_nom[gr=1:N_ext_G])
        @variable(m, LN_s_nom[l=1:N_ext_LN])
        @variable(m, LK_p_nom[l=1:N_ext_LK])
        if investment_type == "continuous"
            @variable(m, LN_inv[l=1:N_ext_LN])
        elseif investment_type == "integer"
            @variable(m, LN_inv[l=1:N_ext_LN], Int)
        elseif investment_type == "binary"
            @variable(m, LN_opt[l=1:N_ext_LN], Bin)
            @variable(m, LN_inv[l=1:N_ext_LN])
        elseif investment_type == "integer_bigm"
            @variable(m, LN_opt[l=1:N_ext_LN,c in candidates[l]], Bin)
        end
        @variable(m, SU_p_nom[s=1:N_ext_SU])
        @variable(m, ST_e_nom[s=1:N_ext_ST])
    end

    if benders != "master"
        @variable(m, G_fix[gr=1:N_fix_G,t=1:Te])
        @variable(m, G_ext[gr=1:N_ext_G,t=1:Te])
        G = [G_fix; G_ext] # G is the concatenated variable array
        @variable(m, LN_fix[l=1:N_fix_LN,t=1:Te])
        @variable(m, LN_ext[l=1:N_ext_LN,t=1:Te])
        LN = [LN_fix; LN_ext]
        @variable(m, LK_fix[l=1:N_fix_LK,t=1:Te])
        @variable(m, LK_ext[l=1:N_ext_LK,t=1:Te])
        LK = [LK_fix; LK_ext]
        @variable(m, SU_dispatch_fix[s=1:N_fix_SU,t=1:Te])
        @variable(m, SU_dispatch_ext[s=1:N_ext_SU,t=1:Te])
        SU_dispatch = [SU_dispatch_fix; SU_dispatch_ext]
        @variable(m, SU_store_fix[s=1:N_fix_SU,t=1:Te])
        @variable(m, SU_store_ext[s=1:N_ext_SU,t=1:Te])
        SU_store = [SU_store_fix; SU_store_ext]
        @variable(m, SU_soc_fix[s=1:N_fix_SU,t=1:Te])
        @variable(m, SU_soc_ext[s=1:N_ext_SU,t=1:Te])
        SU_soc = [SU_soc_fix; SU_soc_ext]
        @variable(m, SU_spill_fix[s=1:N_fix_SU,t=1:Te])
        @variable(m, SU_spill_ext[s=1:N_ext_SU,t=1:Te])
        SU_spill = [SU_spill_fix; SU_spill_ext]
        @variable(m, ST_dispatch_fix[s=1:N_fix_ST,t=1:Te])
        @variable(m, ST_dispatch_ext[s=1:N_ext_ST,t=1:Te])
        ST_dispatch = [ST_dispatch_fix; ST_dispatch_ext]
        @variable(m, ST_store_fix[s=1:N_fix_ST,t=1:Te])
        @variable(m, ST_store_ext[s=1:N_ext_ST,t=1:Te])
        ST_store = [ST_store_fix; ST_store_ext]
        @variable(m, ST_soc_fix[s=1:N_fix_ST,t=1:Te])
        @variable(m, ST_soc_ext[s=1:N_ext_ST,t=1:Te])
        ST_soc = [ST_soc_fix; ST_soc_ext]
        @variable(m, ST_spill_fix[l=1:N_fix_ST,t=1:Te])
        @variable(m, ST_spill_ext[l=1:N_ext_ST,t=1:Te])
        ST_spill = [ST_spill_fix, ST_spill_ext]
        contains(formulation, "angles") ? @variable(m, THETA[1:N,1:Te]) : nothing
    end

    if benders == "master"
        @variable(m, ALPHA[g=1:N_groups]>=0)
    end

    count = 1
    counter = count

    # go through loop only once if full model is built
    if !blockstructure && sn==0
        count = nt
        counter = Trange
        Trange = [Trange]
    elseif !blockstructure && sn>0
        count = nt
        counter = Trange
    end

    for tcurr=Trange

        println("Start building model for snapshot $tcurr.")

    # 1. add all generators to the model
        #println("Adding generators to the model.")

        # 1.1 set different generator types
        generators = network.generators
        fix_gens_b = ((.!generators[:p_nom_extendable]) .& (.!generators[:commitable]))
        ext_gens_b = convert(BitArray, generators[:p_nom_extendable])
        com_gens_b = convert(BitArray, generators[:commitable])

        p_max_pu = get_switchable_as_dense(network, "generators", "p_max_pu")
        p_min_pu = get_switchable_as_dense(network, "generators", "p_min_pu")

        filter_timedependent_extremes!(p_max_pu, 0.01)
        filter_timedependent_extremes!(p_min_pu, 0.01)

        # 1.3a add non-extendable generator variables to the model
        p_min_pu = select_time_dep(network, "generators", "p_min_pu",components=fix_gens_b)
        p_max_pu = select_time_dep(network, "generators", "p_max_pu",components=fix_gens_b)
        p_nom = network.generators[fix_gens_b,:p_nom]

        rf = rescaling_dict[:bounds_G]
        
        if benders != "master"
            if blockstructure || sn>0
                @constraints(m, begin 

                    upper_bounds_G_fix[gr=1:N_fix_G, t=tcurr],
                        rf * G_fix[gr,count] <= rf * p_nom[gr] * p_max_pu(t,gr) 
                
                    lower_bounds_G_fix[gr=1:N_fix_G, t=tcurr],
                        rf * G_fix[gr,count] >= rf * p_nom[gr] * p_min_pu(t,gr) 
                
                end)
            else
                @constraints(m, begin 
                
                    upper_bounds_G_fix[gr=1:N_fix_G, t=tcurr],
                        rf * G_fix[gr,t] <= rf * p_max_pu(t,gr) * p_nom[gr] 
                    
                    lower_bounds_G_fix[gr=1:N_fix_G, t=tcurr],
                        rf * G_fix[gr,t] >= rf * p_nom[gr] * p_min_pu(t,gr)  
                
                end)
            end
        end

        # 1.3b add extendable generator variables to the model
        p_min_pu = select_time_dep(network, "generators", "p_min_pu",components=ext_gens_b)
        p_max_pu = select_time_dep(network, "generators", "p_max_pu",components=ext_gens_b)
        p_nom_min = network.generators[ext_gens_b,:p_nom_min]
        p_nom_max = network.generators[ext_gens_b,:p_nom_max]

        if benders != "slave" && count==nt
            @constraints(m, begin 
                
                lower_bounds_G_p_nom[gr=1:N_ext_G],
                    G_p_nom[gr] >= p_nom_min[gr]
                
                upper_bounds_G_p_nom[gr=1:N_ext_G],
                    G_p_nom[gr] <= p_nom_max[gr]

            end)
        end

        if benders!="master" && benders!="slave"

            @constraints m begin

                lower_bounds_G_ext[t=tcurr,gr=1:N_ext_G],
                    rf * G_ext[gr,t] >= rf * p_min_pu(t,gr) * G_p_nom[gr]

                upper_bounds_G_ext[t=tcurr,gr=1:N_ext_G],
                    rf * G_ext[gr,t] <= rf * p_max_pu(t,gr) * G_p_nom[gr]

            end

        elseif benders == "slave"

            if blockstructure || sn>0
                @constraints(m, begin
                
                    lower_bounds_G_ext[t=tcurr,gr=1:N_ext_G],
                        rf * G_ext[gr,count] >= rf * p_min_pu(t,gr) * generators[ext_gens_b,:p_nom][gr]
                    
                    upper_bounds_G_ext[t=tcurr,gr=1:N_ext_G],
                        rf * G_ext[gr,count] <= rf * p_max_pu(t,gr) * generators[ext_gens_b,:p_nom][gr]

                end)
            else
                @constraints(m, begin

                    lower_bounds_G_ext[t=tcurr,gr=1:N_ext_G],
                        rf*G_ext[gr,t] >= rf*p_min_pu(t,gr)*generators[ext_gens_b,:p_nom][gr]
                
                    upper_bounds_G_ext[t=tcurr,gr=1:N_ext_G],
                        rf*G_ext[gr,t] <= rf*p_max_pu(t,gr)*generators[ext_gens_b,:p_nom][gr]

                end)
            end
        end

        generators = [generators[fix_gens_b,:]; generators[ext_gens_b,:] ] # sort generators the same
        fix_gens_b = ((.!generators[:p_nom_extendable]) .& (.!generators[:commitable]))
        ext_gens_b = convert(BitArray, generators[:p_nom_extendable])
        com_gens_b = convert(BitArray, generators[:commitable])

    # --------------------------------------------------------------------------------------------------------

    # 2. add all lines to the model
        #println("Adding lines to the model.")

        # 2.1 set different lines types
        lines = network.lines
        fix_lines_b = (.!lines[:s_nom_extendable])
        ext_lines_b = .!fix_lines_b

        rf = rescaling_dict[:bounds_LN]

        # 2.3 add variables
        if benders != "master"
            @constraints(m, begin 

                lower_bounds_LN_fix[l=1:N_fix_LN,t=tcurr],
                    rf * LN_fix[l,t] >= rf * (-1) * lines[fix_lines_b,:s_max_pu][l] * lines[fix_lines_b,:s_nom][l]  
            
                upper_bounds_LN_fix[l=1:N_fix_LN,t=tcurr],
                    rf * LN_fix[l,t] <= rf * lines[fix_lines_b,:s_max_pu][l] * lines[fix_lines_b,:s_nom][l]
          
            end)
        end
        
        if benders != "slave" && count==nt
            @constraints(m, begin        
        
                lower_bounds_LN_s_nom[l=1:N_ext_LN],
                    LN_s_nom[l] >= lines[ext_lines_b,:s_nom_min][l] 
                
                upper_bounds_LN_s_nom[l=1:N_ext_LN],
                    LN_s_nom[l] <= lines[ext_lines_b,:s_nom_max][l]
            
            end)
        end
        
        # 2.4 add line constraint for extendable lines
        #println("-- 2.4 add line constraint for extendable lines")

        if benders != "master" && benders != "slave"

            @constraints(m, begin

                upper_bounds_L_ext[t=tcurr,l=1:N_ext_LN], 
                    rf * LN_ext[l,t] <=  rf * lines[ext_lines_b,:s_max_pu][l] * LN_s_nom[l]

                lower_bounds_L_ext[t=tcurr,l=1:N_ext_LN], 
                    rf * LN_ext[l,t] >= rf * (-1) * lines[ext_lines_b,:s_max_pu][l] * LN_s_nom[l]

            end)

        elseif benders == "slave"

            if blockstructure || sn>0
                @constraints(m, begin
    
                    upper_bounds_L_ext[t=tcurr,l=1:N_ext_LN], 
                        rf * LN_ext[l,count] <=  
                        rf * lines[:s_max_pu][l] * lines[ext_lines_b,:s_nom][l]
    
                    lower_bounds_L_ext[t=tcurr,l=1:N_ext_LN],
                        rf * LN_ext[l,count] >=
                        rf * (-1) * lines[:s_max_pu][l] * lines[ext_lines_b,:s_nom][l]
    
                end)
            else
                @constraints(m, begin
    
                    upper_bounds_L_ext[t=tcurr,l=1:N_ext_LN], 
                        rf * LN_ext[l,t] <=  
                        rf * lines[:s_max_pu][l] * lines[ext_lines_b,:s_nom][l]
    
                    lower_bounds_L_ext[t=tcurr,l=1:N_ext_LN],
                        rf * LN_ext[l,t] >=
                        rf * (-1) * lines[:s_max_pu][l] * lines[ext_lines_b,:s_nom][l]
    
                end)
            end

        end

        # 2.5 add integer variables if applicable
        if benders != "slave"  && count==nt

            if investment_type == "continuous"

                @constraint(m, continuous[l=1:N_ext_LN],
                    LN_s_nom[l] ==
                    ( 1.0 + LN_inv[l] / lines[ext_lines_b,:num_parallel][l] ) * 
                    lines[ext_lines_b,:s_nom][l]
                )
                
            elseif investment_type == "integer"

                @constraint(m, integer[l=1:N_ext_LN],
                    LN_s_nom[l] ==
                    ( 1.0 + LN_inv[l] / lines[ext_lines_b,:num_parallel][l] ) *
                    lines[ext_lines_b,:s_nom][l]
                )
                
            elseif investment_type == "binary"

                bigM_default = 1e4
                bigM = min.(lines[ext_lines_b,:s_nom_max],bigM_default)

                @constraints(m, begin

                    binary1[l=1:N_ext_LN],
                        - bigM[l] * ( 1 - LN_opt[l] ) + lines[ext_lines_b,:s_nom_ext_min][l] <= LN_inv[l]

                    binary2[l=1:N_ext_LN],
                        0 <= LN_inv[l] # no reduction of capacity allowed

                    binary3[l=1:N_ext_LN],
                        bigM[l] * LN_opt[l] >= LN_inv[l]

                    binary4[l=1:N_ext_LN],
                        LN_s_nom[l] == 
                        ( 1.0 + LN_inv[l] / lines[ext_lines_b,:num_parallel][l] ) * lines[ext_lines_b,:s_nom][l]

                end)

            elseif investment_type == "integer_bigm"
                
                @constraint(m, integer_bigm_logic[l=1:N_ext_LN], 
                    sum( LN_opt[l,c] for c in candidates[l] ) == 1.0
                )

                @constraint(m, integer_bigm[l=1:N_ext_LN],
                    LN_s_nom[l] ==
                    ( 1 + sum( c * LN_opt[l,c] for c in candidates[l] ) / lines[ext_lines_b,:num_parallel][l] ) *
                    lines[ext_lines_b,:s_nom][l]
                )
                
            end
        end

        lines = [lines[fix_lines_b,:]; lines[ext_lines_b,:]]
        fix_lines_b = (.!lines[:s_nom_extendable])
        ext_lines_b = .!fix_lines_b

    # --------------------------------------------------------------------------------------------------------

    # 3. add all links to the model
        #println("Adding links to the model.")

        rf = rescaling_dict[:bounds_LK]

        # 3.1 set different link types
        links = network.links
        fix_links_b = .!links[:p_nom_extendable]
        ext_links_b = .!fix_links_b

        #  3.3 set link variables
        if benders != "master"

            if blockstructure || sn > 0
                @constraints(m, begin 

                    lower_bounds_LK_fix[l=1:N_fix_LK,t=tcurr],
                        rf * LK_fix[l,count] >= rf * links[fix_links_b, :p_min_pu][l] * links[fix_links_b, :p_nom][l] 
                
                    upper_bounds_LK_fix[l=1:N_fix_LK,t=tcurr],
                        rf * LK_fix[l,count] <= rf * links[fix_links_b, :p_max_pu][l] * links[fix_links_b, :p_nom][l]
            
                end)
            else
                @constraints(m, begin 

                    lower_bounds_LK_fix[l=1:N_fix_LK,t=tcurr],
                        rf * LK_fix[l,t] >= rf * links[fix_links_b, :p_min_pu][l] * links[fix_links_b, :p_nom][l]
                
                    upper_bounds_LK_fix[l=1:N_fix_LK,t=tcurr],
                        rf * LK_fix[l,t] <= rf * links[fix_links_b, :p_max_pu][l] * links[fix_links_b, :p_nom][l]
                
                end)
            end
        end

        if benders != "slave" && count==nt

            @constraints(m, begin 

                lower_bounds_LK_p_nom[l=1:N_ext_LK],
                    LK_p_nom[l] >= links[ext_links_b, :p_nom_min][l]
            
                upper_bounds_LK_p_nom[l=1:N_ext_LK],
                    LK_p_nom[l] <= links[ext_links_b, :p_nom_max][l]
            
            end)
        end

        # # 3.4 set constraints for extendable links
        if benders != "master" && benders != "slave"

            @constraints(m, begin

                lower_bounds_LK_ext[t=tcurr,l=1:N_ext_LK],
                    rf * LK_ext[l,t] >= rf * links[ext_links_b, :p_min_pu][l] * LK_p_nom[l]

                upper_bounds_LK_ext[t=tcurr,l=1:N_ext_LK],
                    rf * LK_ext[l,t] <= rf * links[ext_links_b, :p_max_pu][l] * LK_p_nom[l]
            end)

        elseif benders == "slave"

            if blockstructure || sn > 0

                @constraint(m, lower_bounds_LK_ext[t=tcurr,l=1:N_ext_LK],
                    rf * LK_ext[l,count] >= rf * links[ext_links_b, :p_min_pu][l] * links[ext_links_b,:p_nom][l]
                )

                @constraint(m, upper_bounds_LK_ext[t=tcurr,l=1:N_ext_LK],
                    rf * LK_ext[l,count] <= rf * links[ext_links_b, :p_max_pu][l] * links[ext_links_b,:p_nom][l]
                )

            else 

                @constraint(m, lower_bounds_LK_ext[t=tcurr,l=1:N_ext_LK],
                    rf * LK_ext[l,t] >= rf * links[ext_links_b, :p_min_pu][l] * links[ext_links_b,:p_nom][l]
                )

                @constraint(m, upper_bounds_LK_ext[t=tcurr,l=1:N_ext_LK],
                    rf * LK_ext[l,t] <= rf * links[ext_links_b, :p_max_pu][l] * links[ext_links_b,:p_nom][l]
                )

            end

        end

        links = [links[fix_links_b,:]; links[ext_links_b,:]]
        fix_links_b = .!links[:p_nom_extendable]
        ext_links_b = .!fix_links_b

    # --------------------------------------------------------------------------------------------------------

    # TODO: not made for benders and individual snapshots yet!

    # 4. define storage_units
        #println("Adding storage units to the model.")

        # 4.1 set different storage_units types
        storage_units = network.storage_units
        fix_sus_b = .!storage_units[:p_nom_extendable]
        ext_sus_b = .!fix_sus_b

        inflow = get_switchable_as_dense(network, "storage_units", "inflow")

        #  4.3 set variables
        if benders != "master"

            @constraints(m, begin

                lower_bounds_SU_dispatch_fix[s=1:N_fix_SU,t=counter], 
                    (0 <=  SU_dispatch_fix[s,t])

                upper_bounds_SU_dispatch_fix[s=1:N_fix_SU,t=counter], 
                    (SU_dispatch_fix[s,t] <=
                    (storage_units[fix_sus_b, :p_nom].*storage_units[fix_sus_b, :p_max_pu])[s])
        
                bounds_SU_dispatch_ext[s=1:N_fix_SU,t=counter],
                    SU_dispatch_ext[s,t] >= 0
        
                lower_bounds_SU_store_fix[s=1:N_fix_SU,t=counter],
                    (0 <=  SU_store_fix[s,t])

                upper_bounds_SU_store_fix[s=1:N_fix_SU,t=counter],
                    (SU_store_fix[s,t] <=
                    - (storage_units[fix_sus_b, :p_nom].*storage_units[fix_sus_b, :p_min_pu])[s])
        
                bounds_SU_store_ext[s=1:N_fix_SU,t=counter],
                    SU_store_ext[s,t] >= 0
        
                lower_bounds_SU_soc_fix[s=1:N_fix_SU,t=counter],
                    0 <= SU_soc_fix[s,t]

                upper_bounds_SU_soc_fix[s=1:N_fix_SU,t=counter],
                    SU_soc_fix[s,t] <= 
                    (storage_units[fix_sus_b,:max_hours] .* storage_units[fix_sus_b,:p_nom])[s]
                
                bounds_SU_soc_ext[s=1:N_fix_SU,t=counter],
                    SU_soc_ext[s,t] >= 0

            end)
            
            
            if blockstructure || sn>0
                @constraints(m, begin 
                    lower_bounds_SU_spill_fix[s=1:N_fix_SU,t=tcurr],
                        0 <=  SU_spill_fix[s,count]
                    
                    lower_bounds_SU_spill_ext[s=1:N_fix_SU,t=tcurr],
                        0 <=  SU_spill_ext[s,count]

                    upper_bounds_SU_spill_fix[s=1:N_fix_SU,t=tcurr],
                        SU_spill_fix[s,count] <= inflow[:,fix_sus_b][t,s]
                    
                    upper_bounds_SU_spill_ext[s=1:N_fix_SU,t=tcurr],
                        SU_spill_ext[s,count] <= inflow[:,ext_sus_b][t,s]
                end)
            else
                @constraints(m, begin 
                    lower_bounds_SU_spill_fix[s=1:N_fix_SU,t=tcurr],
                        0 <=  SU_spill_fix[s,t]
                    
                    lower_bounds_SU_spill_ext[s=1:N_fix_SU,t=tcurr],
                        0 <=  SU_spill_ext[s,t]

                    upper_bounds_SU_spill_fix[s=1:N_fix_SU,t=tcurr],
                        SU_spill_fix[s,t] <= inflow[:,fix_sus_b][t,s]
                    
                    upper_bounds_SU_spill_ext[s=1:N_fix_SU,t=tcurr],
                        SU_spill_ext[s,t] <= inflow[:,ext_sus_b][t,s]
                end)
            end
        end
        
        if benders != "slave" && count==nt
            @constraint(m, bounds_SU_p_nom[s=1:N_ext_SU], 
                SU_p_nom[s=1:N_ext_SU] >= 0
            )
        end

        # 4.4 set constraints for extendable storage_units
        if benders != "master" && benders != "slave"
            @constraints(m, begin
                    su_dispatch_limit[t=counter,s=1:N_ext_SU], SU_dispatch_ext[s,t] <= SU_p_nom[s].*storage_units[ext_sus_b, :p_max_pu][s]
                    su_storage_limit[t=counter,s=1:N_ext_SU], SU_store_ext[s,t] <= - SU_p_nom[s].*storage_units[ext_sus_b, :p_min_pu][s]
                    su_capacity_limit[t=counter,s=1:N_ext_SU], SU_soc_ext[s,t] <= SU_p_nom[s].*storage_units[ext_sus_b, :max_hours][s]
            end)
        elseif benders == "slave"
            @constraints(m, begin
                    su_dispatch_limit[t=counter,s=1:N_ext_SU], SU_dispatch_ext[s,t] <= storage_units[ext_sus_b, :p_nom][s].*storage_units[ext_sus_b, :p_max_pu][s]
                    su_storage_limit[t=counter,s=1:N_ext_SU], SU_store_ext[s,t] <= - storage_units[ext_sus_b, :p_nom][s].*storage_units[ext_sus_b, :p_min_pu][s]
                    su_capacity_limit[t=counter,s=1:N_ext_SU], SU_soc_ext[s,t] <= storage_units[ext_sus_b, :p_nom][s].*storage_units[ext_sus_b, :max_hours][s]
            end)
        end 

        storage_units = [storage_units[fix_sus_b,:]; storage_units[ext_sus_b,:]]
        inflow = [inflow[:,fix_sus_b] inflow[:,ext_sus_b]]

        ext_sus_b = BitArray(storage_units[:p_nom_extendable])

        is_cyclic_i = collect(1:N_SU)[BitArray(storage_units[:cyclic_state_of_charge])]
        not_cyclic_i = collect(1:N_SU)[.!storage_units[:cyclic_state_of_charge]]

        if benders != "master"

            @constraints(m, begin

                    su_logic1[s=is_cyclic_i],
                        SU_soc[s,1] == 
                        (SU_soc[s,T]
                        + storage_units[s,:efficiency_store] * SU_store[s,1]
                        - (1./storage_units[s,:efficiency_dispatch]) * SU_dispatch[s,1]
                        + inflow[1,s] - SU_spill[s,1] )

                    su_logic2[s=not_cyclic_i],
                        SU_soc[s,1] == 
                        (storage_units[s,:state_of_charge_initial]
                        + storage_units[s,:efficiency_store] * SU_store[s,1]
                        - (1./storage_units[s,:efficiency_dispatch]) * SU_dispatch[s,1]
                        + inflow[1,s] - SU_spill[s,1])
        
                    su_logic3[s=1:N_SU,t=2:T],
                        SU_soc[s,t] == 
                        (SU_soc[s,t-1]
                        + storage_units[s,:efficiency_store] * SU_store[s,t]
                        - (1./storage_units[s,:efficiency_dispatch]) * SU_dispatch[s,t]
                        + inflow[t,s] - SU_spill[s,t] )
            end)
        end

    # --------------------------------------------------------------------------------------------------------

    # 5. define stores
        #println("Adding stores to the model.")

        # 5.1 set different stores types
        stores = network.stores
        fix_stores_b = .!stores[:e_nom_extendable]
        ext_stores_b = .!fix_stores_b

        inflow = get_switchable_as_dense(network, "stores", "inflow")

        #  5.3 set variables
        if benders != "master"

            @constraints(m, begin

                lower_bounds_ST_dispatch_fix[s=1:N_fix_ST,t=counter],
                    0 <=  ST_dispatch_fix[s,t]

                upper_bounds_ST_dispatch_fix[s=1:N_fix_ST,t=counter],
                    (ST_dispatch_fix[s,t] <=
                    (stores[fix_stores_b, :e_nom].*stores[fix_stores_b, :e_max_pu])[s])

                bounds_ST_dispatch_ext[s=1:N_fix_ST,t=counter],
                    ST_dispatch_ext[s,t] >= 0

                lower_bounds_ST_store_fix[s=1:N_fix_ST,t=counter],
                    0 <=  ST_store_fix[s,t]

                upper_bounds_ST_store_fix[s=1:N_fix_ST,t=counter],
                    (ST_store_fix[s,t] <=
                    - (stores[fix_stores_b, :e_nom].*stores[fix_stores_b, :e_min_pu])[s])

                bounds_ST_store_ext[s=1:N_fix_ST,t=counter],
                    ST_store_ext[s,t] >= 0

                lower_bounds_ST_soc_fix[s=1:N_fix_ST,t=counter],
                    0 <= ST_soc_fix[s,t]

                upper_bounds_ST_soc_fix[s=1:N_fix_ST,t=counter],
                    ST_soc_fix[s,t] <= 
                    (stores[fix_stores_b,:max_hours] .* stores[fix_stores_b,:e_nom])[s]

                bounds_ST_soc_ext[s=1:N_fix_ST,t=counter],
                    ST_soc_ext[s,t] >= 0

            end)
            
            
            if blockstructure || sn>0
                @constraints(m, begin 
                    lower_bounds_ST_spill_fix[s=1:N_fix_ST,t=tcurr],
                        0 <=  ST_spill_fix[s,count]
    
                    lower_bounds_ST_spill_ext[s=1:N_fix_ST,t=tcurr],
                        0 <=  ST_spill_ext[s,count]

                    upper_bounds_ST_spill_fix[s=1:N_fix_ST,t=tcurr],
                        ST_spill_fix[s,count] <= 
                        inflow[:,fix_stores_b][s,t]
    
                    upper_bounds_ST_spill_ext[s=1:N_fix_ST,t=tcurr],
                        ST_spill_ext[s,count] <= 
                        inflow[:,ext_stores_b][s,t]
                end)
            else 
                @constraints(m, begin 
                    lower_bounds_ST_spill_fix[s=1:N_fix_ST,t=tcurr],
                        0 <=  ST_spill_fix[s,t]
    
                    lower_bounds_ST_spill_ext[s=1:N_fix_ST,t=tcurr],
                        0 <=  ST_spill_ext[s,t]
                  
                    upper_bounds_ST_spill_fix[s=1:N_fix_ST,t=tcurr],
                        ST_spill_fix[s,t] <= 
                        inflow[:,fix_stores_b][s,t]
    
                    upper_bounds_ST_spill_ext[s=1:N_fix_ST,t=tcurr],
                        ST_spill_ext[s,t] <= 
                        inflow[:,ext_stores_b][s,t]
                end)
            end
        end
        
        if benders != "slave" && count==nt
            @constraint(m, bounds_ST_e_nom[s=1:N_ext_ST],
                ST_e_nom[s] >= 0
            )
        end

        # 5.4 set constraints for extendable stores
        if benders != "master" && benders != "slave"
            @constraints(m, begin
                    st_dispatch_limit[t=counter,s=1:N_ext_ST], 
                        ST_dispatch_ext[s,t] <= ST_e_nom[s].*stores[ext_stores_b, :e_max_pu][s]
                    st_storage_limit[t=counter,s=1:N_ext_ST], 
                        ST_store_ext[s,t] <= - ST_e_nom[s].*stores[ext_stores_b, :e_min_pu][s]
                    st_capacity_limit[t=counter,s=1:N_ext_ST], 
                        ST_soc_ext[s,t] <= ST_e_nom[s].*stores[ext_stores_b, :max_hours][s]
            end)
        elseif benders == "slave"
            @constraints(m, begin
                    st_dispatch_limit[t=counter,s=1:N_ext_ST],
                        ST_dispatch_ext[s,t] <= 
                        stores[ext_stores_b, :e_nom][s].*stores[ext_stores_b, :e_max_pu][s]
                    st_storage_limit[t=counter,s=1:N_ext_ST],
                        ST_store_ext[s,t] <= 
                        - stores[ext_stores_b, :e_nom][s].*stores[ext_stores_b, :e_min_pu][s]
                    st_capacity_limit[t=counter,s=1:N_ext_ST],
                        ST_soc_ext[s,t] <= 
                        stores[ext_stores_b, :e_nom][s].*stores[ext_stores_b, :max_hours][s]
            end)
        end

        # 5.5 set charging constraint
        stores = [stores[fix_stores_b,:]; stores[ext_stores_b,:]]
        inflow = [inflow[:,fix_stores_b] inflow[:,ext_stores_b]]

        ext_stores_b = BitArray(stores[:e_nom_extendable])

        is_cyclic_i = collect(1:N_ST)[BitArray(stores[:cyclic_state_of_charge])]
        not_cyclic_i = collect(1:N_ST)[.!stores[:cyclic_state_of_charge]]

        if benders != "master"

            @constraints(m, begin

                    st_logic1[s=is_cyclic_i,t=1], 
                        ST_soc[s,t] == 
                        (ST_soc[s,T]
                        + stores[s,:efficiency_store] * ST_store[s,t]
                        - stores[s,:efficiency_dispatch] * ST_dispatch[s,t]
                        + inflow[t,s] - ST_spill[s,t])
        
                    st_logic2[s=not_cyclic_i,t=1], 
                        ST_soc[s,t] == 
                        (stores[s,:state_of_charge_initial]
                        + stores[s,:efficiency_store] * ST_store[s,t]
                        - stores[s,:efficiency_dispatch] * ST_dispatch[s,t]
                        + inflow[t,s] - ST_spill[s,t])
        
                    st_logic3[s=is_cyclic_i,t=2:T], 
                        ST_soc[s,t] == 
                        (ST_soc[s,t-1]
                        + stores[s,:efficiency_store] * ST_store[s,t]
                        - stores[s,:efficiency_dispatch] * ST_dispatch[s,t]
                        + inflow[t,s] - ST_spill[s,t])

            end)
        end

    # --------------------------------------------------------------------------------------------------------

    # 6. power flow formulations

        #println("Adding power flow formulation $formulation to the model.")

        rf = rescaling_dict[:flows]

        if benders != "master"

            # a.1 linear angles formulation
            if formulation == "angles_linear"

                # load data in correct order
                loads = network.loads_t["p"][:,Symbol.(network.loads[:name])]

                if sn>0
                    @constraint(m, nodal[t=counter,n=1:N],
    
                            sum(G[findin(generators[:bus], [reverse_busidx[n]]), count])
                            + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                            .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,count])
                            + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,count])
    
                            - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                            - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,count])
                            - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,count])
    
                            == sum(LN[findin(lines[:bus0], [reverse_busidx[n]]),count])
                            - sum(LN[findin(lines[:bus1], [reverse_busidx[n]]),count]) 
                    )

                else

                    @constraint(m, nodal[t=counter,n=1:N],
    
                            sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                            + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                            .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                            + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])
    
                            - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                            - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                            - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])
    
                            == sum(LN[findin(lines[:bus0], [reverse_busidx[n]]),t])
                            - sum(LN[findin(lines[:bus1], [reverse_busidx[n]]),t]) 
                    )
                end

                sn>0 ? cnt=count : cnt=counter

                @constraint(m, flows[t=cnt,l=1:L], 
                    rf * LN[l, t] == 
                    rf * lines[:x_pu][l]^(-1) *     
                    ( THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t] )
                )

                @constraint(m, slack[t=cnt], THETA[1,t] == 0 )

            elseif formulation == "angles_linear_integer_bigm"

                bigm_upper = bigm(:flows_upper, network)
                bigm_lower = bigm(:flows_lower, network)

                # load data in correct order
                loads = network.loads_t["p"][:,Symbol.(network.loads[:name])]

                if benders!="master" && benders!="slave"

                    @constraint(m, nodal[t=counter,n=1:N], (

                        sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                        + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                            .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                        + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                        - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                        - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                        - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                        == sum(LN[findin(lines[:bus0], [reverse_busidx[n]]),t])
                        - sum(LN[findin(lines[:bus1], [reverse_busidx[n]]),t]) ))

                    @constraint(m, flows_upper[l=(sum(fix_lines_b)+1):(sum(ext_lines_b)+sum(fix_lines_b)),c in candidates[l],t=counter],
                        rf * (
                            ( 1 + c / lines[:num_parallel][l] ) * lines[:x_pu][l]^(-1) * 
                            ( THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t] )
                            - LN[l,t]
                        ) >= rf * ( LN_opt[l,c] - 1 ) * bigm_upper
                    )
                    
                    @constraint(m, flows_lower[l=(sum(fix_lines_b)+1):(sum(ext_lines_b)+sum(fix_lines_b)),c in candidates[l],t=counter],
                        rf * (
                            ( 1 + c / lines[:num_parallel][l] ) * lines[:x_pu][l]^(-1) * 
                            ( THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t] )
                            - LN[l,t]
                        ) >= rf * ( 1 - LN_opt[l,c] ) * bigm_lower    
                    )

                    @constraint(m, flows_fix[l=1:(sum(fix_lines_b)), t=counter],

                        rf * LN[l, t] == 
                        rf * lines[:x_pu][l]^(-1) * ( THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t] ) 
                    
                    )   

                    @constraint(m, slack[t=counter], THETA[1,t] == 0)

                elseif benders=="slave"

                    if blockstructure || sn>0

                        @constraint(m, nodal[t=counter,n=1:N], (

                            sum(G[findin(generators[:bus], [reverse_busidx[n]]), count])
                            + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                                .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,count])
                            + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,count])

                            - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                            - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,count])
                            - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,count])

                            == sum(LN[findin(lines[:bus0], [reverse_busidx[n]]),count])
                            - sum(LN[findin(lines[:bus1], [reverse_busidx[n]]),count]) ))

                        @constraint(m, flows_upper[l=(sum(fix_lines_b)+1):(sum(ext_lines_b)+sum(fix_lines_b)),c in candidates[l],t=counter],
                            rf * (
                                ( 1 + c / lines[:num_parallel][l] ) * lines[:x_pu][l]^(-1) * 
                                ( THETA[busidx[lines[:bus0][l]], count] - THETA[busidx[lines[:bus1][l]], count] ) 
                                - LN[l,count]
                            ) >= rf * bigm_upper * (c == 0 ? 0 : -1)
                        )
        
                        @constraint(m, flows_lower[l=(sum(fix_lines_b)+1):(sum(ext_lines_b)+sum(fix_lines_b)),c in candidates[l],t=counter],
                            rf * (
                                ( 1 + c / lines[:num_parallel][l] ) * lines[:x_pu][l]^(-1) ) * 
                                ( THETA[busidx[lines[:bus0][l]], count] - THETA[busidx[lines[:bus1][l]], count] ) 
                                - LN[l,count]
                            ) <= rf * bigm_lower * (c == 0 ? 0 : 1)
                        )

                        @constraint(m, flows_nonext[l=1:(sum(fix_lines_b)), t=counter],
                            LN[l, count] == 
                            lines[:x_pu][l]^(-1) * ( THETA[busidx[lines[:bus0][l]], count] - THETA[busidx[lines[:bus1][l]], count] ) 
                        ) 
                        
                        @constraint(m, slack[t=counter], THETA[1,count] == 0)

                    else

                        @constraint(m, nodal[t=counter,n=1:N], (

                            sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                            + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                                .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                            + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                            - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                            - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                            - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                            == sum(LN[findin(lines[:bus0], [reverse_busidx[n]]),t])
                            - sum(LN[findin(lines[:bus1], [reverse_busidx[n]]),t]) ))

                        @constraint(m, flows_upper[l=(sum(fix_lines_b)+1):(sum(ext_lines_b)+sum(fix_lines_b)),c in candidates[l],t=tcurr],
                            rf * (
                                ( 1 + c / lines[:num_parallel][l] ) * lines[:x_pu][l]^(-1) * 
                                ( THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t] )
                                - LN[l,t]
                            ) >= rf * bigm_upper*(c == 0 ? 0 : -1)
                        )
        
                        @constraint(m, flows_lower[l=(sum(fix_lines_b)+1):(sum(ext_lines_b)+sum(fix_lines_b)),c in candidates[l],t=tcurr],
                            rf * (
                                ( 1 + c / lines[:num_parallel][l] ) *lines[:x_pu][l]^(-1) *
                                ( THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t] )
                                - LN[l,t]
                            ) <= rf * bigm_lower * (c == 0 ? 0 : 1)
                        )

                        @constraint(m, flows_nonext[l=1:(sum(fix_lines_b)), t=tcurr],
                            LN[l, t] == lines[:x_pu][l]^(-1) *
                            ( THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t] ) 
                        )  

                        @constraint(m, slack[t=tcurr], THETA[1,t] == 0)

                    end
                end   

            elseif formulation == "angles_bilinear"

                # cannot be solved by Gurobi, needs Ipopt!
                # needs investment_type defined!
                # requires starting reactance and line capacity!

                # load data in correct order
                loads = network.loads_t["p"][:,Symbol.(network.loads[:name])]

                @constraint(m, nodal[t=counter,n=1:N], (

                        sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                        + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                                .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                        + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                        - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                        - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                        - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                        == sum(LN[findin(lines[:bus0], [reverse_busidx[n]]),t])
                        - sum(LN[findin(lines[:bus1], [reverse_busidx[n]]),t]) )        
                )

                @NLconstraint(m, flows_ext[t=counter,l=(sum(fix_lines_b)+1):(sum(ext_lines_b)+sum(fix_lines_b))],
                    LN[l, t] ==  
                    (1 + LN_inv[l] / lines[:num_parallel][l]) * lines[:x_pu][l]^(-1) *
                    (THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t] ) 
                )

                @constraint(m, flows_nonext[t=counter,l=1:(sum(fix_lines_b))],
                    LN[l, t] == lines[:x_pu][l]^(-1) *
                    ( THETA[busidx[lines[:bus0][l]], t] - THETA[busidx[lines[:bus1][l]], t] ) 
                )                                                   

                @constraint(m, slack[t=counter], THETA[1,t] == 0)

            # b.1 linear kirchhoff formulation
            elseif formulation == "kirchhoff_linear"

                #load data in correct order
                loads = network.loads_t["p"][:,Symbol.(network.loads[:name])]

                @constraint(m, balance[t=counter,n=1:N], (
                    sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                    + sum(LN[ findin(lines[:bus1], [reverse_busidx[n]]) ,t])
                    + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                        .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                    + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                    - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                    - sum(LN[ findin(lines[:bus0], [reverse_busidx[n]]) ,t])
                    - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                    - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                    == 0 )
                )

                # Might be nessecary to loop over all subgraphs as
                # for (sn, sub) in enumerate(weakly_connected_components(g))
                #     g_sub = induced_subgraph(g, sub)[1]

                append_idx_col!(lines)
                (branches, var, attribute) = (lines, LN, :x)
                cycles = get_cycles(network)
                if ndims(cycles)<2
                    cycles = [cycle for cycle in cycles if length(cycle)>2]
                else
                    cycles = [cycles[i,:] for i in 1:size(cycles)[1]]
                end
                if length(cycles)>0
                    cycles_branch = Array{Int64,1}[]
                    directions = Array{Float64,1}[]
                    for cyc=1:length(cycles)
                        push!(cycles_branch,Int64[])
                        push!(directions,Float64[])
                        for bus=1:length(cycles[cyc])
                            bus0 = cycles[cyc][bus]
                            bus1 = cycles[cyc][(bus)%length(cycles[cyc])+1]
                            try
                                push!(cycles_branch[cyc],branches[((branches[:bus0].==reverse_busidx[bus0])
                                            .&(branches[:bus1].==reverse_busidx[bus1])),:idx][1] )
                                push!(directions[cyc], 1.)
                            catch y
                                if isa(y, BoundsError)
                                    push!(cycles_branch[cyc], branches[((branches[:bus0].==reverse_busidx[bus1])
                                                    .&(branches[:bus1].==reverse_busidx[bus0])),:idx][1] )
                                    push!(directions[cyc], -1.)
                                else
                                    return y
                                end
                            end
                        end
                    end
                    if attribute==:x
                        @constraint(m, line_cycle_constraint[t=counter,c=1:length(cycles_branch)] ,
                                dot(directions[c] .* lines[cycles_branch[c], :x_pu],
                                    LN[cycles_branch[c],t]) == 0
                        )
                    # elseif attribute==:r
                    #     @constraint(m, link_cycle_constraint[c=1:length(cycles_branch), t=counter] ,
                    #             dot(directions[c] .* links[cycles_branch[c], :r]/380. , LK[cycles_branch[c],t]) == 0)
                    end
                end

            # b.2 bilinear kirchhoff formulation (steps derived from original s_nom and x)
            elseif formulation == "kirchhoff_bilinear"

                # cannot be solved by Gurobi, needs Ipopt!
                # needs investment_type defined!
                # requires starting reactance and line capacity!

                #load data in correct order
                loads = network.loads_t["p"][:,Symbol.(network.loads[:name])]

                @constraint(m, balance[t=counter,n=1:N], (

                    sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                    + sum(LN[ findin(lines[:bus1], [reverse_busidx[n]]) ,t])
                    + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                        .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                    + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                    - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                    - sum(LN[ findin(lines[:bus0], [reverse_busidx[n]]) ,t])
                    - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                    - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                    == 0 )
                )

                # Might be nessecary to loop over all subgraphs as
                # for (sn, sub) in enumerate(weakly_connected_components(g))
                #     g_sub = induced_subgraph(g, sub)[1]

                append_idx_col!(lines)
                (branches, var, attribute) = (lines, LN, :x)
                cycles = get_cycles(network)
                if ndims(cycles)<2
                    cycles = [cycle for cycle in cycles if length(cycle)>2]
                else
                    cycles = [cycles[i,:] for i in 1:size(cycles)[1]]
                end
                if length(cycles)>0
                    cycles_branch = Array{Int64,1}[]
                    directions = Array{Float64,1}[]
                    for cyc=1:length(cycles)
                        push!(cycles_branch,Int64[])
                        push!(directions,Float64[])
                        for bus=1:length(cycles[cyc])
                            bus0 = cycles[cyc][bus]
                            bus1 = cycles[cyc][(bus)%length(cycles[cyc])+1]
                            try
                                push!(cycles_branch[cyc],branches[((branches[:bus0].==reverse_busidx[bus0])
                                            .&(branches[:bus1].==reverse_busidx[bus1])),:idx][1] )
                                push!(directions[cyc], 1.)
                            catch y
                                if isa(y, BoundsError)
                                    push!(cycles_branch[cyc], branches[((branches[:bus0].==reverse_busidx[bus1])
                                                    .&(branches[:bus1].==reverse_busidx[bus0])),:idx][1] )
                                    push!(directions[cyc], -1.)
                                else
                                    return y
                                end
                            end
                        end
                    end
                    if attribute==:x
                        @NLconstraint(m, line_cycle_constraint[t=counter,c=1:length(cycles_branch)] ,
                                sum(      directions[c][l] 
                                        * lines[cycles_branch[c], :x_pu][l] 
                                        * (1+LN_inv[cycles_branch[c]][l] / lines[cycles_branch[c], :num_parallel][l])^(-1)
                                        #--* (1+LN_inv[cycles_branch[c]][l])^(-1)
                                        * LN[cycles_branch[c],t][l]
                                for l=1:length(directions[c])
                                )               
                                    == 0)
                    # elseif attribute==:r
                    #     @constraint(m, link_cycle_constraint[c=1:length(cycles_branch), t=counter] ,
                    #             dot(directions[c] .* links[cycles_branch[c], :r]/380. , LK[cycles_branch[c],t]) == 0)
                    end
                end

            # c.1 linear ptdf formulation
            elseif formulation == "ptdf"

                ptdf = ptdf_matrix(network)

                #load data in correct order
                loads = network.loads_t["p"][:,Symbol.(network.loads[:name])]

                @constraint(m, flows[t=counter,l=1:L],
                    sum( ptdf[l,n]
                    * (   sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                        + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                            .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                        + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                        - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                        - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                        - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])
                        )
                    for n in 1:N
                    ) == LN[l,t] 
                )

                @constraint(m, balance[t=counter],
                    sum( (   sum(G[findin(generators[:bus], [reverse_busidx[n]]), t])
                        + sum(links[findin(links[:bus1], [reverse_busidx[n]]),:efficiency]
                            .* LK[ findin(links[:bus1], [reverse_busidx[n]]) ,t])
                        + sum(SU_dispatch[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])

                        - row_sum(loads[t,findin(network.loads[:bus],[reverse_busidx[n]])],1)
                        - sum(LK[ findin(links[:bus0], [reverse_busidx[n]]) ,t])
                        - sum(SU_store[ findin(storage_units[:bus], [reverse_busidx[n]]) ,t])
                        ) for n in 1:N
                    ) == 0 
                )

            else
                #println("The formulation $formulation is not implemented.")
            end
        end
    

    # --------------------------------------------------------------------------------------------------------

    # 7. set global_constraints

        # 7.1 carbon constraint
        if count==nt

            #println("Adding global CO2 constraints to the model.")
    
            carrier_index(carrier) = findin(generators[:carrier], carrier)
    
            if benders != "master" && sn==0
                if nrow(network.global_constraints)>0 && in("co2_limit", network.global_constraints[:name])
                    co2_limit = network.global_constraints[network.global_constraints[:name].=="co2_limit", :constant]
                    println("CO2_limit is $(co2_limit[1]) t")
                    nonnull_carriers = network.carriers[network.carriers[:co2_emissions].!=0, :][:name]
                    @constraint(m, co2limit, sum(sum(network.snapshots[:weightings][t]*dot(1./generators[carrier_index(nonnull_carriers) , :efficiency],
                                G[carrier_index(nonnull_carriers),t]) for t=1:T)
                                * select_names(network.carriers, [carrier])[:co2_emissions]
                                for carrier in network.carriers[:name]) .<=  co2_limit)
                end
            end
    
            # 7.2 limit on transmission expansion volume
            # sum of capacity expansion times length over all lines <= limit in MWkm unit
            if benders != "slave"
                if nrow(network.global_constraints)>0 && in(true, network.lines[:s_nom_extendable]) && in("mwkm_limit", network.global_constraints[:name])
                    mwkm_limit = network.global_constraints[network.global_constraints[:name].=="mwkm_limit", :constant]
                    #println("Line expansion limit is $(mwkm_limit[1]) times current MWkm")
                    @constraint(m, mwkmlimit, 
                        dot(LN_s_nom,lines[:length]) <= mwkm_limit[1] * dot(lines[ext_lines_b,:s_nom],lines[:length])
                    )
                end
            end
    
            # TODO: fix
            # 7.3 specified percentage of renewable energy generation
            # sum of renewable generation =/>= percentage * sum of total load
            if benders != "master" && sn==0
                if nrow(network.global_constraints)>0 && in("restarget", network.global_constraints[:name])
                    restarget = network.global_constraints[network.global_constraints[:name].=="restarget", :constant]
                    #println("Target share of renewable energy is $(restarget[1]*100) %")
                    null_carriers = network.carriers[network.carriers[:co2_emissions].==0,:][:name]
                    @constraint(m, restarget,
                        sum(sum(network.snapshots[:weightings][t]*G[carrier_index(null_carriers),t] for t=1:T))
                        .>= restarget * sum(network.snapshots[:weightings][t]*sum(convert(Array,network.loads_t["p_set"][t,2:end])) for t=1:T)
                    )
                end
            end
    
            # 7.4 fake specified percentages of renewable energy generation (if there is no/little curtailment)
            if benders != "slave"
                if nrow(network.global_constraints)>0 && in("approx_restarget", network.global_constraints[:name])
    
                    N_loads = size(network.loads_t["p_set"])[2]
                    approx_restarget = network.global_constraints[network.global_constraints[:name].=="approx_restarget", :constant]
    
                    #println("Target share of renewable energy is $(approx_restarget[1]*100) %")
    
                    # needed to take out biomass since it is not time dependent
                    null_carriers = network.carriers[(network.carriers[:co2_emissions].==0) .& (network.carriers[:name] .!= "biomass"),:][:name]
                    
                    ren_gens_b = [in(i,carrier_index(null_carriers)) ? true : false for i=1:size(generators)[1]]
                    fix_ren_gens_b = .!generators[:p_nom_extendable][ren_gens_b]
                    ext_ren_gens_b = .!fix_ren_gens_b
                    
                    ren_gens_b_orig = [in(i,findin(network.generators[:carrier], null_carriers)) ? true : false for i=1:size(network.generators)[1]]
                    fix_ren_gens_b_orig = .!network.generators[:p_nom_extendable][ren_gens_b_orig]
                    ext_ren_gens_b_orig = .!fix_ren_gens_b_orig
                    
                    fix_ren_gens = generators[fix_gens_b .& ren_gens_b,:]
                    ext_ren_gens = generators[ext_gens_b .& ren_gens_b,:]
                    
                    def_p_max_pu_ext = 8760.0*ext_ren_gens[:p_max_pu]
                    def_p_max_pu_fix = 8760.0*fix_ren_gens[:p_max_pu]
                    
                    exist_fix_ren_gens = maximum(fix_ren_gens_b_orig)

                    p_max_pu_full = network.generators_t["p_max_pu"][:,2:end]
                    exist_fix_ren_gens ? p_max_pu_fix = convert(Array,p_max_pu_full[:,fix_ren_gens_b_orig]) : nothing
                    p_max_pu_ext = convert(Array,p_max_pu_full[:,ext_ren_gens_b_orig])
    
                    loc_fix = findin(fix_ren_gens[:name], string.(p_max_pu_full.colindex.names))
                    loc_ext = findin(ext_ren_gens[:name], string.(p_max_pu_full.colindex.names))
    
                    loc_fix_b = [in(i,loc_fix) ? true : false for i=1:length(def_p_max_pu_fix)]
                    loc_ext_b = [in(i,loc_ext) ? true : false for i=1:length(def_p_max_pu_ext)]
                    
                    exist_fix_ren_gens ? sum_of_p_max_pu_fix = sum(network.snapshots[:weightings][t]*p_max_pu_fix[t,:] for t=1:T)  : nothing
                    sum_of_p_max_pu_ext = sum(network.snapshots[:weightings][t]*p_max_pu_ext[t,:] for t=1:T)
    
                    exist_fix_ren_gens ? def_p_max_pu_fix[loc_fix_b] .= sum_of_p_max_pu_fix : nothing
                    def_p_max_pu_ext[loc_ext_b] .= sum_of_p_max_pu_ext

                    rf = rescaling_dict[:approx_restarget]

                    if exist_fix_ren_gens

                        @constraint(m, approx_restarget,
                            rf * (dot(def_p_max_pu_fix,fix_ren_gens[:p_nom]) + dot(def_p_max_pu_ext,G_p_nom))
                            >= rf * approx_restarget[1] * sum(network.snapshots[:weightings][t]*network.loads_t["p_set"][t,n] for t=1:T for n=2:N_loads)
                        )
                        
                    else
                        @constraint(m, approx_restarget,
                            rf * dot(def_p_max_pu_ext,G_p_nom)
                            >= rf * approx_restarget[1] * sum(network.snapshots[:weightings][t]*network.loads_t["p_set"][t,n] for t=1:T for n=2:N_loads)
                        )
                    end
                end
            end
        end
        

    # --------------------------------------------------------------------------------------------------------

    # 8. set objective function

        if count==nt
            #println("Adding objective to the model.")

            if benders!="master" && benders!="slave"
    
                @objective(m, Min,
                        sum(network.snapshots[:weightings][t] * dot(generators[:marginal_cost], G[:,t]) for t=1:Te)
                        + dot(generators[ext_gens_b,:capital_cost], G_p_nom[:] )
                        + dot(generators[fix_gens_b,:capital_cost], generators[fix_gens_b,:p_nom])
    
                        + dot(lines[ext_lines_b,:capital_cost], LN_s_nom[:])
                        + dot(lines[fix_lines_b,:capital_cost], lines[fix_lines_b,:s_nom])
    
                        + dot(links[ext_links_b,:capital_cost], LK_p_nom[:])
                        + dot(links[fix_links_b,:capital_cost], links[fix_links_b,:p_nom])
    
                        + sum(network.snapshots[:weightings][t] * dot(storage_units[:marginal_cost], SU_dispatch[:,t]) for t=1:Te)
                        + dot(storage_units[ext_sus_b, :capital_cost], SU_p_nom[:])
                        + dot(storage_units[fix_sus_b,:capital_cost], storage_units[fix_sus_b,:p_nom])
    
                        + sum(network.snapshots[:weightings][t] * dot(stores[:marginal_cost], ST_dispatch[:,t]) for t=1:Te)
                        + dot(stores[ext_stores_b, :capital_cost], ST_e_nom[:])
                        + dot(stores[fix_stores_b,:capital_cost], stores[fix_stores_b,:e_nom])
                )
    
            elseif benders == "master"
    
                @objective(m, Min,
                    dot(generators[ext_gens_b,:capital_cost], G_p_nom[:] )
                    + dot(generators[fix_gens_b,:capital_cost], generators[fix_gens_b,:p_nom])
            
                    + dot(lines[ext_lines_b,:capital_cost], LN_s_nom[:])
                    + dot(lines[fix_lines_b,:capital_cost], lines[fix_lines_b,:s_nom])
            
                    + dot(links[ext_links_b,:capital_cost], LK_p_nom[:])
                    + dot(links[fix_links_b,:capital_cost], links[fix_links_b,:p_nom])
            
                    + dot(storage_units[ext_sus_b, :capital_cost], SU_p_nom[:])
                    + dot(storage_units[fix_sus_b,:capital_cost], storage_units[fix_sus_b,:p_nom])
            
                    + dot(stores[ext_stores_b, :capital_cost], ST_e_nom[:])
                    + dot(stores[fix_stores_b,:capital_cost], stores[fix_stores_b,:e_nom])
            
                    + sum( ALPHA[g] for g=1:N_groups )
                )

            elseif benders == "slave"
                
                if sn>0
                    @objective(m, Min,
                              sum( network.snapshots[:weightings][t] * dot(generators[:marginal_cost], G[:,count]) for t=counter )
                            + sum( network.snapshots[:weightings][t] * dot(storage_units[:marginal_cost], SU_dispatch[:,count]) for t=counter )
                            + sum( network.snapshots[:weightings][t] * dot(stores[:marginal_cost], ST_dispatch[:,count]) for t=counter )
                    )
                else
                    @objective(m, Min,
                              sum( network.snapshots[:weightings][t] * dot(generators[:marginal_cost], G[:,t]) for t=counter )
                            + sum( network.snapshots[:weightings][t] * dot(storage_units[:marginal_cost], SU_dispatch[:,t]) for t=counter )
                            + sum( network.snapshots[:weightings][t] * dot(stores[:marginal_cost], ST_dispatch[:,t]) for t=counter )
                    )
                end
    
            else
                error()
            end
        end

        count += 1
        counter = count

    end # for loop

    return m

end
