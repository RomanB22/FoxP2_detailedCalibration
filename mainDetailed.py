import bluepyopt as bpopt
import bluepyopt.ephys as ephys
import matplotlib.pyplot as plt

plotMorpho=False
plotResponses=True
offspringSize=2
maxGenerations=2

morphoFile = './morphology/detailedMorpho.swc' # detailedMorpho threeCompartmental
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
    'g_pas.all': 7.5e-05,
    'e_pas.all': -80,
    'cm.all': 1,
    'cm.basal': 2,
    'cm.apical': 2,
    # Sodium currents
    'gbar_NaTs.somatic': 1,
    'gbar_NaTs.basal': 2,
    'gbar_NaTs.apical': 2,
    'gbar_NaTs.axonal': 10,
    'gbar_Nap.somatic': 1e-7,
    'gbar_Nap.basal': 1e-7,
    'gbar_Nap.apical': 1e-7,
    'gbar_Nap.axonal': 1e-7,
    # Potassium currents
    'gbar_Kv3_1.somatic': 3e-1,
    'gbar_Kv3_1.basal': 3e-1,
    'gbar_Kv3_1.apical': 3e-1,
    'gbar_Kv3_1.axonal': 3e-1,
    'gbar_K_T.somatic': 1e-3,
    'gbar_K_T.basal': 1e-3,
    'gbar_K_T.apical': 1e-3,
    'gbar_K_T.axonal': 1e-3,
    'gbar_K_P.somatic': 1e-2,
    'gbar_K_P.basal': 1e-2,
    'gbar_K_P.apical': 1e-2,
    'gbar_K_P.axonal': 1e-2,
    'gbar_SK.somatic': 3e-3,
    'gbar_SK.basal': 3e-3,
    'gbar_SK.apical': 3e-3,
    'gbar_SK.axonal': 3e-3,
    # I-HCN current
    'gbar_Ih.somatic': 1e-4,
    'gbar_Ih.basal': 1e-4,
    'gbar_Ih.apical': 1e-4,
    'gbar_Ih.axonal': 1e-4,
    # I-M current
    'gbar_Im.somatic': 1e-3,
    'gbar_Im.basal': 1e-3,
    'gbar_Im.apical': 1e-3,
    'gbar_Im.axonal': 1e-3,
    # CaHVA current
    'gbar_Ca_HVA.somatic': 1e-4,
    'gbar_Ca_HVA.basal': 1e-4,
    'gbar_Ca_HVA.apical': 1e-4,
    'gbar_Ca_HVA.axonal': 1e-4,
    # CaLVA current
    'gbar_Ca_LVA.somatic': 1e-3,
    'gbar_Ca_LVA.basal': 1e-3,
    'gbar_Ca_LVA.apical': 1e-3,
    'gbar_Ca_LVA.axonal': 1e-3
}

def plot_responses(responses, filename='./figures/optResults/responses.png'):
    fig, axes = plt.subplots(len(responses), figsize=(10,8))
    for index, (resp_name, response) in enumerate(responses.items()):
        axes[index].plot(response['time'], response['voltage'], label=resp_name)
        axes[index].set_title(resp_name)
    fig.tight_layout()
    fig.savefig(filename)
if plotResponses: 
    release_responses = evaluator.run_protocols(protocols=fitness_protocols.values(), param_values=release_params)
    plot_responses(release_responses, filename='./figures/optResults/Original_responses.png')
    plt.show()

# opt = bpopt.optimisations.DEAPOptimisation(
#     evaluator=evaluator,
#     offspring_size=offspringSize)

opt = bpopt.deapext.optimisationsCMA.DEAPOptimisationCMA(
    evaluator=evaluator,
    offspring_size=offspringSize,
    selector_name="multi_objective") #single_objective multi_objective

final_pop, halloffame, log, hist = opt.run(max_ngen=maxGenerations, cp_filename='checkpoints/checkpoint.pkl')

for i in range(len(halloffame)):
    best_params = evaluator.param_dict(halloffame[i])
    print('Best %d solution: \n' % i + str(best_params))

    best_responses = evaluator.run_protocols(protocols=fitness_protocols.values(), param_values=best_params)
    if plotResponses: plot_responses(best_responses, filename='./figures/optResults/Best_responses_%i.png' % i)
