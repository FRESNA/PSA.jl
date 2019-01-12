using JuMP
using MathProgBase

# TODO: neglects storage units and storages at the moment

"""Solves a linear optimal power flow with Benders Decomposition using lazy constraints."""
function run_lazybenders_lopf(network, solver; 
    formulation::String = "angles_linear",
    investment_type::String = "continuous",
    rescaling::Bool = false,
    split_subproblems::Bool = false,
    individualcuts::Bool = false,
    tolerance::Float64 = 100.0,
    mip_gap::Float64 = 1e-8,
    bigM::Float64 = 1e12,
    update_x::Bool = false
    )

    #################
    # sanity checks #
    #################

    investment_type=="integer_bigm" && formulation!="angles_linear_integer_bigm" ? error("Investment type $investment_type requires formulation <angles_linear_integer_bigm>, not $formulation!") : nothing
    update_x && investment_type=="integer_bigm" ? error("Integer Big-M formulation does not support update_x option") : nothing

    ###################
    # precalculations #
    ###################

    calculate_dependent_values!(network)
    T = nrow(network.snapshots)
    ext_gens_b = ext_gens_b = convert(BitArray, network.generators[:p_nom_extendable])
    ext_lines_b = convert(BitArray, network.lines[:s_nom_extendable])
    ext_links_b = convert(BitArray, network.links[:p_nom_extendable])
    N_ext_G = sum(ext_gens_b)
    N_ext_LN = sum(ext_lines_b)
    N_ext_LK = sum(ext_links_b)
    individualcuts ? N_cuts = T : N_cuts = 1

    p_min_pu = select_time_dep(network, "generators", "p_min_pu",components=ext_gens_b)
    p_max_pu = select_time_dep(network, "generators", "p_max_pu",components=ext_gens_b)
    filter_timedependent_extremes!(p_max_pu, 0.01)
    filter_timedependent_extremes!(p_min_pu, 0.01)

    x_orig = network.lines[:x]
    x_pu_orig = network.lines[:x_pu]

    bigm_upper = bigm(:flows_upper, network)
    bigm_lower = bigm(:flows_lower, network)

    rf = 1
    rf_dict = rescaling_factors(rescaling)

    coupled_slave_constrs = [
        :lower_bounds_G_ext,
        :upper_bounds_G_ext,
        :lower_bounds_LN_ext,
        :upper_bounds_LN_ext,
        ]

    if N_ext_LK > 0
        push!(coupled_slave_constrs, :lower_bounds_LK_ext)
        push!(coupled_slave_constrs, :upper_bounds_LK_ext)
    end

    if investment_type=="integer_bigm"
        push!(coupled_slave_constrs, :flows_upper)
        push!(coupled_slave_constrs, :flows_lower)
    end

    ###################################
    # build master and slave problems # 
    ###################################

    # build master
    model_master = build_lopf(network, solver,
        benders="master",
        investment_type=investment_type,
        rescaling=rescaling,
        N_cuts=N_cuts
    )

    if !update_x

        # build subproblems
        models_slave = JuMP.Model[]
        if !split_subproblems
            push!(models_slave, 
                build_lopf(network, solver, 
                    benders="slave",
                    formulation=formulation,
                    rescaling=rescaling
                )
            )
        else
            for t=1:T
                push!(models_slave, 
                    build_lopf(network, solver, 
                        benders="slave",
                        formulation=formulation,
                        rescaling=rescaling,
                        snapshot_number=t
                    )
                )
            end
        end
        N_slaves = length(models_slave)
    
        uncoupled_slave_constrs = setdiff(getconstraints(first(models_slave)), coupled_slave_constrs)
        
    end
    
    #################################
    # benders cut callback function #
    #################################

    function benderscut(cb)

        # get results of master problem
        objective_master_current = getobjectivevalue(model_master)
        G_p_nom_current = getvalue(model_master[:G_p_nom])
        LN_s_nom_current = getvalue(model_master[:LN_s_nom])
        N_ext_LK > 0 ? LK_p_nom_current = getvalue(model_master[:LK_p_nom]) : nothing
        investment_type=="integer_bigm" ? LN_opt_current = getvalue(model_master[:LN_opt]) : LN_inv_current = getvalue(model_master[:LN_inv])

        ALPHA_current = sum(getvalue(model_master[:ALPHA][g]) for g=1:N_cuts)

        candidates = line_extensions_candidates(network)

        if update_x

            # update line reactances (pu) according to master values
            for l=1:nrow(network.lines)
                if network.lines[:s_nom_extendable][l]
                    network.lines[:x_pu][l] = x_pu_orig[l] * ( 1 + LN_inv_current[l] / network.lines[:num_parallel][l] )
                    network.lines[:x][l] = x_orig[l] * ( 1 + LN_inv_current[l] / network.lines[:num_parallel][l] )
                end
            end
            
            # build subproblems
            models_slave = JuMP.Model[]
            if !split_subproblems
                push!(models_slave, 
                    build_lopf(network, solver, 
                        benders="slave",
                        formulation=formulation,
                    )
                )
            else
                for t=1:T
                    push!(models_slave, 
                        build_lopf(network, solver, 
                            benders="slave",
                            formulation=formulation,
                            snapshot_number=t
                        )
                    )
                end
            end
            N_slaves = length(models_slave)
        
            uncoupled_slave_constrs = setdiff(getconstraints(first(models_slave)), coupled_slave_constrs)
        end


        # adapt RHS of slave models with solution of master model

        model_slave = models_slave[1]
        
        for t=1:T

            split_subproblems ? model_slave = models_slave[t] : nothing


            # RHS generators

            rf = rf_dict[:bounds_G]

            for gr=1:N_ext_G

                JuMP.setRHS(
                    model_slave[:lower_bounds_G_ext][t,gr],
                    rf * G_p_nom_current[gr] * p_min_pu(t,gr)
                )

                JuMP.setRHS(
                    model_slave[:upper_bounds_G_ext][t,gr],
                    rf * G_p_nom_current[gr] * p_max_pu(t,gr)
                )

            end

            # RHS lines

            rf = rf_dict[:bounds_LN]

            for l=1:N_ext_LN

                JuMP.setRHS(
                    model_slave[:lower_bounds_LN_ext][t,l],
                    rf * (-1) * network.lines[ext_lines_b,:s_max_pu][l] * LN_s_nom_current[l]
                )
                
                JuMP.setRHS(
                    model_slave[:upper_bounds_LN_ext][t,l],
                    rf * network.lines[ext_lines_b,:s_max_pu][l] * LN_s_nom_current[l]
                )

            end

            # RHS links

            if N_ext_LK> 0

                rf = rf_dict[:bounds_LK]
    
                for l=1:N_ext_LK
    
                    JuMP.setRHS(
                        model_slave[:lower_bounds_LK_ext][t,l],
                        rf * LK_p_nom_current[l] * network.links[ext_links_b, :p_min_pu][l]
                    )
                    
                    JuMP.setRHS(
                        model_slave[:upper_bounds_LK_ext][t,l],
                        rf * LK_p_nom_current[l] * network.links[ext_links_b, :p_max_pu][l]
                    )
    
                end

            end

            # RHS flows
                    
            rf = rf_dict[:flows]

            for l=1:N_ext_LN
                
                if investment_type=="integer_bigm"
                   
                    for c in candidates[l]

                        rhs = rf * ( LN_opt_current[l,c] - 1 ) * bigm_upper[l]
                        if (rhs < 1e-4 && rhs > -1e-4)
                            rhs = 0.0
                        end

                        JuMP.setRHS(model_slave[:flows_upper][l,c,t],rhs)
                        
                        rhs = rf * ( 1 - LN_opt_current[l,c] ) * bigm_lower[l]
                        if (rhs < 1e-4 && rhs > -1e-4)
                            rhs = 0.0
                        end

                        JuMP.setRHS(model_slave[:flows_lower][l,c,t], rhs)

                    end

                end

            end

        end

        # solve slave problems
        statuses_slave = Symbol[]
        for i=1:N_slaves
            push!(statuses_slave, solve(models_slave[i])) 
        end
        status_slave = minimum(statuses_slave)

        # print iis if infeasible
        for i in length(statuses_slave)
            if statuses_slave[i] == :Infeasible && typeof(solver) == Gurobi.GurobiSolver
                println("WARNING: Subproblem $i is infeasible. The IIS is:")
                println(get_iis(models_slave[i]))
            end
        end

        # get results of slave problems
        objective_slave_current = sum(getobjectivevalue(models_slave[i]) for i=1:N_slaves)
        duals_lower_bounds_G_ext = getduals(models_slave, :lower_bounds_G_ext)
        duals_upper_bounds_G_ext = getduals(models_slave, :upper_bounds_G_ext)
        duals_lower_bounds_LN_ext = getduals(models_slave, :lower_bounds_LN_ext)
        duals_upper_bounds_LN_ext = getduals(models_slave, :upper_bounds_LN_ext)

        if N_ext_LK > 0 
            duals_lower_bounds_LK_ext = getduals(models_slave, :lower_bounds_LK_ext)
            duals_upper_bounds_LK_ext = getduals(models_slave, :upper_bounds_LK_ext)
        end
        
        if investment_type=="integer_bigm"
            duals_flows_upper = getduals_flows(models_slave, :flows_upper)
            duals_flows_lower = getduals_flows(models_slave, :flows_lower)
        end

        # go through cases of subproblems
        if (status_slave == :Optimal && 
            abs(objective_slave_current - ALPHA_current) <= tolerance)

            println("Optimal solution of the original problem found")
            println("The optimal objective value is ", objective_master_current)
            return
        end

        if !individualcuts
            T_range = [1:T] 
        else 
            T_range = 1:T
        end

        for T_range_curr = T_range

            N_slaves == 1 ? slave_id = 1 : slave_id = T_range_curr
            individualcuts ? cut_id = T_range_curr : cut_id = 1

            # calculate coupled cut components

            cut_G_lower = sum(  
                    duals_lower_bounds_G_ext[t,gr] * rf_dict[:bounds_G] * p_min_pu(t,gr) * model_master[:G_p_nom][gr]
                for t=T_range_curr for gr=1:N_ext_G)
            
            cut_G_upper = sum(  
                    duals_upper_bounds_G_ext[t,gr] * rf_dict[:bounds_G] * p_max_pu(t,gr) * model_master[:G_p_nom][gr]
                for t=T_range_curr for gr=1:N_ext_G)

            cut_LN_lower = sum(  
                    duals_lower_bounds_LN_ext[t,l] * rf_dict[:bounds_LN] * (-1) * network.lines[ext_lines_b,:s_max_pu][l] * model_master[:LN_s_nom][l]
                for t=T_range_curr for l=1:N_ext_LN)

            cut_LN_upper = sum(  
                    duals_upper_bounds_LN_ext[t,l] * rf_dict[:bounds_LN] * network.lines[ext_lines_b,:s_max_pu][l] * model_master[:LN_s_nom][l]
                for t=T_range_curr for l=1:N_ext_LN)

            cut_LK_lower = 0
            cut_LK_upper = 0

            if N_ext_LK > 0
                cut_LK_lower = sum( 
                        duals_lower_bounds_LK_ext[t,l] * rf_dict[:bounds_LK] * network.links[ext_links_b,:p_min_pu][l] * model_master[:LK_p_nom][l]
                    for t=T_range_curr for l=1:N_ext_LK)
    
                cut_LK_upper = sum( 
                        duals_upper_bounds_LK_ext[t,l] * rf_dict[:bounds_LK] * network.links[ext_links_b,:p_max_pu][l] * model_master[:LK_p_nom][l]
                    for t=T_range_curr for l=1:N_ext_LK)
            end

            cut_G = cut_G_lower + cut_G_upper
            cut_LN = cut_LN_lower + cut_LN_upper
            cut_LK = cut_LK_lower + cut_LK_upper

            if investment_type=="integer_bigm"

                cut_flows_lower = sum(  
                        duals_flows_lower[slave_id][l,c,t] * rf_dict[:flows] * ( 1 - model_master[:LN_opt][l,c] ) * bigm_lower[l] 
                    for t=T_range_curr for l=1:N_ext_LN for c in candidates[l])

                cut_flows_upper = sum(  
                        duals_flows_upper[slave_id][l,c,t] * rf_dict[:flows] * ( model_master[:LN_opt][l,c] - 1 ) * bigm_upper[l] 
                    for t=T_range_curr for l=1:N_ext_LN for c in candidates[l])

                cut_flows = cut_flows_lower + cut_flows_upper

            end

            # calculate uncoupled cut components
            cut_const = get_benderscut_constant(models_slave[slave_id],uncoupled_slave_constrs)

            if (status_slave == :Optimal &&
                abs(objective_slave_current - ALPHA_current) > tolerance)

                rf = rf_dict[:benderscut]

                if investment_type!="integer_bigm"

                    @lazyconstraint(cb, rf * model_master[:ALPHA][cut_id] >= 
                        rf * ( cut_G + cut_LN + cut_LK + cut_const ) 
                    )

                else

                    @lazyconstraint(cb, rf * model_master[:ALPHA][cut_id] >= 
                        rf * ( cut_G + cut_LN + cut_LK + cut_flows + cut_const ) 
                    )

                end
            end
            
            if status_slave == :Infeasible

                if investment_type!="integer_bigm"

                    @lazyconstraint(cb, 0 >= 
                        rf * ( cut_G + cut_LN + cut_LK + cut_const ) 
                    )
                
                else
                    
                    @lazyconstraint(cb, 0 >=
                        rf * ( cut_G + cut_LN + cut_LK + cut_flows + cut_const ) 
                    )

                end
                
            end

        end

    end # function benderscut

    ###################################
    # solve master with lazy callback #
    ###################################

    addlazycallback(model_master, benderscut)
    status_master = solve(model_master);

    # print iis if infeasible
    if status_master == :Infeasible && typeof(solver) == Gurobi.GurobiSolver
        println("ERROR: Master problem is infeasible. The IIS is:")
        println(get_iis(model_master))
    end

    ##################
    # output results #
    ##################

    write_optimalsolution(network, model_master; sm=models_slave, joint=false)
    return (model_master, models_slave);

end