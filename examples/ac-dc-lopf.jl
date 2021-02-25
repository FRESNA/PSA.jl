
# Cf. PyPSA example https://github.com/PyPSA/PyPSA/tree/master/examples/ac-dc-meshed

import PSA

using Clp
using JuMP

network = PSA.import_network("/home/tom/fias/lib/pypsa/examples/ac-dc-meshed/ac-dc-data/")

solver = Clp.Optimizer

m = PSA.lopf(network, solver)

println(objective_value(m))

println(network.generators_t["p"])
