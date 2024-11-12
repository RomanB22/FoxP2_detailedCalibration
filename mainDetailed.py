import bluepyopt as bpopt
import bluepyopt.ephys as ephys
import matplotlib.pyplot as plt

plotMorpho=False
plotResponses=True
offspringSize=100
maxGenerations=200

morphoFile = './morphology/threeCompartmental.swc'
workDir = './config_3Comp/'

if plotMorpho:
    import neurom
    from neurom.view import plot_morph, plot_morph3d
    plot_morph3d(neurom.load_morphology(morphoFile)) #
    plt.show()

    # from neurom import features
    # import numpy as np
    # m = neurom.load_morphology(morphoFile)
    # SomaSurf = features.get('soma_surface_area', m)
    # NeuritesSurf = np.sum(features.get('total_area_per_neurite', m))
    # print(features.get('total_area_per_neurite', m),features.get('total_length_per_neurite', m))
    # print("Soma surface area: " + str(SomaSurf) + " um", '\n', "Neuron surface area: " + str(SomaSurf+NeuritesSurf) + " um")
    quit()
morphology = ephys.morphologies.NrnFileMorphology(morphology_path=morphoFile, do_replace_axon=True)

import json
param_configs = json.load(open(workDir+'parameters.json'))
print([param_config['param_name'] for param_config in param_configs])

import l5pc_model
parameters = l5pc_model.define_parameters()
print('\n'.join('%s' % param for param in parameters))

mechanisms = l5pc_model.define_mechanisms()
print('\n'.join('%s' % mech for mech in mechanisms))

l5pc_cell = ephys.models.CellModel('l5pc', morph=morphology, mechs=mechanisms, params=parameters)
print(l5pc_cell)

param_names = [param.name for param in l5pc_cell.params.values() if not param.frozen]

proto_configs = json.load(open(workDir+'protocols.json'))
print(proto_configs)

import l5pc_evaluator
fitness_protocols = l5pc_evaluator.define_protocols()
print('\n'.join('%s' % protocol for protocol in fitness_protocols.values()))

feature_configs = json.load(open(workDir+'features.json'))
print(feature_configs)

fitness_calculator = l5pc_evaluator.define_fitness_calculator(fitness_protocols)
# print(fitness_calculator)

sim = ephys.simulators.NrnSimulator()

evaluator = ephys.evaluators.CellEvaluator(
        cell_model=l5pc_cell,
        param_names=param_names,
        fitness_protocols=fitness_protocols,
        fitness_calculator=fitness_calculator,
        sim=sim)

release_params = {
    'gIhbar_Ih.apical': 0,
    'gIhbar_Ih.somatic': 0,
    'g_pas.all': 7.5e-05,
    'e_pas.all': -80,
    'cm.all': 1,
    'gNaTs2_tbar_NaTs2_t.apical': 0.026145,
    'gSKv3_1bar_SKv3_1.apical': 0.004226,
    'gImbar_Im.apical': 0.000143,
    'gNaTa_tbar_NaTa_t.axonal': 3.137968,
    'gK_Tstbar_K_Tst.axonal': 0.089259,
    'gamma_CaDynamics_E2.axonal': 0.002910,
    'gNap_Et2bar_Nap_Et2.axonal': 0.006827,
    'gSK_E2bar_SK_E2.axonal': 0.007104,
    'gCa_HVAbar_Ca_HVA.axonal': 0.000990,
    'gK_Pstbar_K_Pst.axonal': 0.973538,
    'gSKv3_1bar_SKv3_1.axonal': 1.021945,
    'decay_CaDynamics_E2.axonal': 287.198731,
    'gCa_LVAstbar_Ca_LVAst.axonal': 0.008752,
    'gamma_CaDynamics_E2.somatic': 0.000609,
    'gSKv3_1bar_SKv3_1.somatic': 0.303472,
    'gSK_E2bar_SK_E2.somatic': 0.008407,
    'gCa_HVAbar_Ca_HVA.somatic': 0.000994,
    'gNaTs2_tbar_NaTs2_t.somatic': 0.983955,
    'decay_CaDynamics_E2.somatic': 210.485284,
    'gCa_LVAstbar_Ca_LVAst.somatic': 0.000333
}

def plot_responses(responses, filename='./figures/responses.png'):
    fig, axes = plt.subplots(len(responses), figsize=(10,8))
    for index, (resp_name, response) in enumerate(sorted(responses.items())):
        axes[index].plot(response['time'], response['voltage'], label=resp_name)
        axes[index].set_title(resp_name)
    fig.tight_layout()
    fig.savefig(filename)
if plotResponses: 
    release_responses = evaluator.run_protocols(protocols=fitness_protocols.values(), param_values=release_params)
    plot_responses(release_responses, filename='./figures/Original_responses.png')

# opt = bpopt.optimisations.DEAPOptimisation(
#     evaluator=evaluator,
#     offspring_size=offspringSize)

opt = bpopt.deapext.optimisationsCMA.DEAPOptimisationCMA(
    evaluator=evaluator,
    offspring_size=offspringSize,
    selector_name="multi_objective") #single_objective multi_objective

final_pop, halloffame, log, hist = opt.run(max_ngen=maxGenerations, cp_filename='checkpoints/checkpoint.pkl')

best_params = evaluator.param_dict(halloffame[0])
print(best_params)

best_responses = evaluator.run_protocols(protocols=fitness_protocols.values(), param_values=best_params)
if plotResponses: plot_responses(best_responses, filename='./figures/Best_responses.png')

